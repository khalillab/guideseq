#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

guideseq.py
===========

V2

serves as the wrapper for all guideseq pipeline

"""

import os
import sys
import yaml
import argparse
import traceback
# Set up logger
import log
logger = log.createCustomLogger('root')

from alignReads import alignReads
from filterBackgroundSites import filterBackgroundSites
from umi import demultiplex, umitag, consolidate
from visualization import visualizeOfftargets
import identifyOfftargetSites
import validation
from tabulate import tabulate

DEFAULT_DEMULTIPLEX_MIN_READS = 10000
DEFAULT_WINDOW_SIZE = 25
DEFAULT_MAX_SCORE = 7

CONSOLIDATE_MIN_QUAL = 15
CONSOLIDATE_MIN_FREQ = 0.9

def get_parameters(manifest_data):
    default_yaml = os.path.dirname(os.path.realpath(__file__)) + "/default.yaml"
    with open(default_yaml, 'r') as f:
        default = yaml.load(f)
    with open(manifest_data, 'r') as f:
        return_dict = yaml.load(f)
    default['analysis_folder'] = os.getcwd()
    validation.validateManifest(return_dict)
    return_dict['parameters'] = {}
    for p in default:
        if not p in return_dict:
            if p == "samples":
                logger.error("No samples are found in the yaml file, please provide samples (fastq) to start with!")
                exit()
            elif p == "undemultiplexed":
                logger.error("No undemultiplexed samples are found in the yaml file, please provide raw fastq to start with!")
                exit()
            else:
                
                return_dict['parameters'][p] = default[p]
    for s in return_dict['samples']:
        try:
            return_dict['samples'][s]['control_barcode1']
        except:
            return_dict['samples'][s]['control_barcode1'] = "AAAAAAAA"
        try:
            return_dict['samples'][s]['control_barcode2']
        except:
            return_dict['samples'][s]['control_barcode2'] = "AAAAAAAA"
    return return_dict

class GuideSeq:

    def __init__(self):
        pass

    def parseManifest(self, manifest_path):
        logger.info('Loading manifest...')

        try:
            manifest_data = get_parameters(manifest_path)
            self.undemultiplexed = manifest_data['undemultiplexed']
            self.samples = manifest_data['samples']
            self.parameters = manifest_data['parameters']
            print(tabulate([(k,v) for k,v in self.parameters.items()])) 
            print(tabulate([(k,v) for k,v in self.undemultiplexed.items()])) 
            print(tabulate([(k,v) for k,v in self.samples[list(self.samples.keys())[0]].items()])) 

        except Exception as e:
            logger.error('Error demultiplexing reads.')
            logger.error(traceback.format_exc())
            sys.exit()

        logger.info('Successfully loaded manifest.')

    def demultiplex(self):

        logger.info('Demultiplexing undemultiplexed files...')
        # Take our two barcodes and concatenate them
        swapped_sample_barcodes = {} # Name is confusing? Yichao notes
        for sample in self.samples:
            barcode1 = self.samples[sample]['barcode1']
            barcode2 = self.samples[sample]['barcode2']
            # if control barcodes do not exist, get_parameter function will assign AAAAAAAA as the barcode
            control_barcode1 = self.samples[sample]['control_barcode1']
            control_barcode2 = self.samples[sample]['control_barcode2']
            barcode = barcode1[1:8] + barcode2[1:8]
            control_barcode = control_barcode1[1:8] + control_barcode2[1:8]
            swapped_sample_barcodes[barcode] = sample
            swapped_sample_barcodes[control_barcode] = "control_"+sample

        try:
            demultiplex.demultiplex(self.undemultiplexed['forward'],
                                    self.undemultiplexed['reverse'],
                                    self.undemultiplexed['index1'],
                                    self.undemultiplexed['index2'],
                                    swapped_sample_barcodes,
                                    os.path.join(self.parameters['analysis_folder'], 'demultiplexed'),
                                    min_reads=self.parameters['demultiplex_min_reads'])

            self.demultiplexed = {}
            for sample in self.samples:
                self.demultiplexed[sample] = {}
                self.demultiplexed[sample]['read1'] = os.path.join(self.parameters['analysis_folder'], 'demultiplexed', sample + '.r1.fastq')
                self.demultiplexed[sample]['read2'] = os.path.join(self.parameters['analysis_folder'], 'demultiplexed', sample + '.r2.fastq')
                self.demultiplexed[sample]['index1'] = os.path.join(self.parameters['analysis_folder'], 'demultiplexed', sample + '.i1.fastq')
                self.demultiplexed[sample]['index2'] = os.path.join(self.parameters['analysis_folder'], 'demultiplexed', sample + '.i2.fastq')
                sample = "control_"+sample
                self.demultiplexed[sample] = {}
                self.demultiplexed[sample]['read1'] = os.path.join(self.parameters['analysis_folder'], 'demultiplexed', sample + '.r1.fastq')
                self.demultiplexed[sample]['read2'] = os.path.join(self.parameters['analysis_folder'], 'demultiplexed', sample + '.r2.fastq')
                self.demultiplexed[sample]['index1'] = os.path.join(self.parameters['analysis_folder'], 'demultiplexed', sample + '.i1.fastq')
                self.demultiplexed[sample]['index2'] = os.path.join(self.parameters['analysis_folder'], 'demultiplexed', sample + '.i2.fastq')

            logger.info('Successfully demultiplexed reads.')
        except Exception as e:
            logger.error('Error demultiplexing reads.')
            logger.error(traceback.format_exc())
            quit()

    def umitag(self):
        logger.info('umitagging reads...')

        try:
            self.umitagged = {}
            for sample in self.samples:
                try:
                    self.umitagged[sample] = {}
                    self.umitagged[sample]['read1'] = os.path.join(self.parameters['analysis_folder'], 'umitagged', sample + '.r1.umitagged.fastq')
                    self.umitagged[sample]['read2'] = os.path.join(self.parameters['analysis_folder'], 'umitagged', sample + '.r2.umitagged.fastq')

                    umitag.umitag(self.demultiplexed[sample]['read1'],
                                  self.demultiplexed[sample]['read2'],
                                  self.demultiplexed[sample]['index1'],
                                  self.demultiplexed[sample]['index2'],
                                  self.umitagged[sample]['read1'],
                                  self.umitagged[sample]['read2'],
                                  os.path.join(self.parameters['analysis_folder'], 'umitagged'))

                    sample = "control_"+sample
                    self.umitagged[sample] = {}
                    self.umitagged[sample]['read1'] = os.path.join(self.parameters['analysis_folder'], 'umitagged', sample + '.r1.umitagged.fastq')
                    self.umitagged[sample]['read2'] = os.path.join(self.parameters['analysis_folder'], 'umitagged', sample + '.r2.umitagged.fastq')

                    umitag.umitag(self.demultiplexed[sample]['read1'],
                                  self.demultiplexed[sample]['read2'],
                                  self.demultiplexed[sample]['index1'],
                                  self.demultiplexed[sample]['index2'],
                                  self.umitagged[sample]['read1'],
                                  self.umitagged[sample]['read2'],
                                  os.path.join(self.parameters['analysis_folder'], 'umitagged'))
                except:
                    logger.error(f'Failed for sample {sample}')

            logger.info('Successfully umitagged reads.')
        except Exception as e:
            logger.error('Error umitagging')
            logger.error(traceback.format_exc())
            quit()

    def consolidate(self, min_freq=CONSOLIDATE_MIN_FREQ, min_qual=CONSOLIDATE_MIN_QUAL):
        logger.info('Consolidating reads...')

        try:
            self.consolidated = {}
            for sample in self.samples:
                try:
                    self.consolidated[sample] = {}
                    self.consolidated[sample]['read1'] = os.path.join(self.parameters['analysis_folder'], 'consolidated', sample + '.r1.consolidated.fastq')
                    self.consolidated[sample]['read2'] = os.path.join(self.parameters['analysis_folder'], 'consolidated', sample + '.r2.consolidated.fastq')

                    consolidate.consolidate(self.umitagged[sample]['read1'], self.consolidated[sample]['read1'], min_qual, min_freq)
                    consolidate.consolidate(self.umitagged[sample]['read2'], self.consolidated[sample]['read2'], min_qual, min_freq)

                    sample = "control_"+sample
                    self.consolidated[sample] = {}
                    self.consolidated[sample]['read1'] = os.path.join(self.parameters['analysis_folder'], 'consolidated', sample + '.r1.consolidated.fastq')
                    self.consolidated[sample]['read2'] = os.path.join(self.parameters['analysis_folder'], 'consolidated', sample + '.r2.consolidated.fastq')

                    consolidate.consolidate(self.umitagged[sample]['read1'], self.consolidated[sample]['read1'], min_qual, min_freq)
                    consolidate.consolidate(self.umitagged[sample]['read2'], self.consolidated[sample]['read2'], min_qual, min_freq)

                except:
                    logger.error(f'Failed for sample {sample}')

            logger.info('Successfully consolidated reads.')
        except Exception as e:
            logger.error('Error umitagging')
            logger.error(traceback.format_exc())
            quit()

    def alignReads(self):
        logger.info('Aligning reads...')

        try:
            self.aligned = {}
            for sample in self.samples:
                try:
                    sample_alignment_path = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.sam')
                    alignReads(self.parameters['bwa'],
                               self.parameters['reference_genome'],
                               self.consolidated[sample]['read1'],
                               self.consolidated[sample]['read2'],
                               sample_alignment_path,njobs=self.parameters['njobs'])
                    self.aligned[sample] = sample_alignment_path
                    sample = "control_"+sample
                    sample_alignment_path = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.sam')
                    alignReads(self.parameters['bwa'],
                               self.parameters['reference_genome'],
                               self.consolidated[sample]['read1'],
                               self.consolidated[sample]['read2'],
                               sample_alignment_path,njobs=self.parameters['njobs'])
                    self.aligned[sample] = sample_alignment_path
                except:
                    logger.error(f'Failed for sample {sample}')    
                logger.info('Finished aligning reads to genome.')

        except Exception as e:
            logger.error('Error aligning')
            logger.error(traceback.format_exc())
            quit()

    def identifyOfftargetSites(self):
        logger.info('Identifying offtarget sites...')

        try:
            self.identified = {}

            # Identify offtarget sites for each sample
            for sample in self.samples:
                try:
                    # Prepare sample annotations
                    sample_data = self.samples[sample]
                    annotations = {}
                    annotations['Description'] = sample_data['description']
                    annotations['Targetsite'] = sample
                    print ("Using control primer",sample_data['control_primer'])


                    samfile = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.sam')

                    self.identified[sample] = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.txt')

                    identifyOfftargetSites.analyze(samfile, self.parameters['reference_genome'], self.identified[sample], annotations,
                                                   self.parameters['window_size'], self.parameters['max_score'], sample_data['control_primer'])

                    sample = "control_"+sample

                    print ("Using control primer",sample_data['control_primer'])


                    samfile = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.sam')

                    self.identified[sample] = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.txt')

                    identifyOfftargetSites.analyze(samfile, self.parameters['reference_genome'], self.identified[sample], annotations,
                                                   self.parameters['window_size'], self.parameters['max_score'], sample_data['control_primer'])
                except:
                    logger.error(f'Failed for sample {sample}')   

            logger.info('Finished identifying offtarget sites.')

        except Exception as e:
            logger.error('Error identifying offtarget sites.')
            logger.error(traceback.format_exc())
            quit()

    def visualize(self):
        logger.info('Visualizing off-target sites')


        for sample in self.samples: ## 3/6/2020 Yichao solved: visualization stopped when one sample failed
            try:
                infile = self.identified[sample]
                outfile = os.path.join(self.parameters['analysis_folder'], 'visualization', sample + '_offtargets')
                try:
                    self.PAM
                    visualizeOfftargets(infile, outfile, title=sample,PAM=self.PAM)
                except:
                    visualizeOfftargets(infile, outfile, title=sample,PAM="NGG")
            except Exception as e:
                logger.error('Error visualizing off-target sites: %s'%(sample))
                logger.error(traceback.format_exc())
            try:
                sample = "control_"+sample
                infile = self.identified[sample]
                outfile = os.path.join(self.parameters['analysis_folder'], 'visualization', sample + '_offtargets')
                try:
                    self.PAM
                    visualizeOfftargets(infile, outfile, title=sample,PAM=self.PAM)
                except:
                    visualizeOfftargets(infile, outfile, title=sample,PAM="NGG")
            except Exception as e:
                logger.error('Error visualizing off-target sites: %s'%(sample))
                logger.error(traceback.format_exc())
        logger.info('Finished visualizing off-target sites')

def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')

    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    all_parser.add_argument('--identifyAndFilter', action='store_true', default=False)
    all_parser.add_argument('--skip_demultiplex', action='store_true', default=False)
    all_parser.add_argument('--test', action='store_true', default=False)

    parallel_parser = subparsers.add_parser('parallel', help='Run all steps of the pipeline in parallel')
    parallel_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -R rusage[mem=60000] -P GUIDEV2 -q priority')
    # parallel_parser.add_argument('--run', '-r', help='Specify which steps of pipepline to run (all, align, identify, visualize, variants)', default='all')


    demultiplex_parser = subparsers.add_parser('demultiplex', help='Demultiplex undemultiplexed FASTQ files')
    demultiplex_parser.add_argument('--manifest', '-m', help='Specify the manifest path', required=True)

    umitag_parser = subparsers.add_parser('umitag', help='UMI tag demultiplexed FASTQ files for consolidation')
    umitag_parser.add_argument('--read1', required=True)
    umitag_parser.add_argument('--read2', required=True)
    umitag_parser.add_argument('--index1', required=True)
    umitag_parser.add_argument('--index2', required=True)
    umitag_parser.add_argument('--outfolder', required=True)

    consolidate_parser = subparsers.add_parser('consolidate', help='Consolidate UMI tagged FASTQs')
    consolidate_parser.add_argument('--read1', required=True)
    consolidate_parser.add_argument('--read2', required=True)
    consolidate_parser.add_argument('--outfolder', required=True)
    consolidate_parser.add_argument('--min_quality', required=False, type=float)
    consolidate_parser.add_argument('--min_frequency', required=False, type=float)

    align_parser = subparsers.add_parser('align', help='Paired end read mapping to genome')
    align_parser.add_argument('--bwa', required=True)
    align_parser.add_argument('--genome', required=True)
    align_parser.add_argument('--read1', required=True)
    align_parser.add_argument('--read2', required=True)
    align_parser.add_argument('--outfolder', required=True)

    identify_parser = subparsers.add_parser('identify', help='Identify GUIDE-seq offtargets')
    identify_parser.add_argument('--aligned', required=True)
    identify_parser.add_argument('--genome', required=True)
    identify_parser.add_argument('--outfolder', required=True)
    identify_parser.add_argument('--target_sequence', required=True)
    identify_parser.add_argument('--description', required=False)
    identify_parser.add_argument('--max_score', required=False, type=int, default=7)
    identify_parser.add_argument('--window_size', required=False, type=int, default=25)

    filter_parser = subparsers.add_parser('filter', help='Filter identified sites from control sites')
    filter_parser.add_argument('--bedtools', required=True)
    filter_parser.add_argument('--identified', required=True)
    filter_parser.add_argument('--background', required=True)
    filter_parser.add_argument('--outfolder', required=True)

    visualize_parser = subparsers.add_parser('visualize', help='Visualize off-target sites')
    visualize_parser.add_argument('--infile', required=True)
    visualize_parser.add_argument('--outfolder', required=True)
    visualize_parser.add_argument('--title', required=False)

    return parser.parse_args()


def main():
    args = parse_args()

    if args.command == 'all':
        if args.test:
            g = GuideSeq()
            g.parseManifest(args.manifest)
            exit()
        if args.identifyAndFilter:
            try:
                g = GuideSeq()
                g.parseManifest(args.manifest)
                g.identifyOfftargetSites()
                g.visualize()

            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        elif args.skip_demultiplex:
            try:
                g = GuideSeq()
                g.parseManifest(args.manifest)
                g.demultiplexed = {}
                for sample in g.samples:
                    g.demultiplexed[sample] = {}
                    g.demultiplexed[sample]['read1'] = os.path.join(g.analysis_folder, 'demultiplexed', sample + '.r1.fastq')
                    g.demultiplexed[sample]['read2'] = os.path.join(g.analysis_folder, 'demultiplexed', sample + '.r2.fastq')
                    g.demultiplexed[sample]['index1'] = os.path.join(g.analysis_folder, 'demultiplexed', sample + '.i1.fastq')
                    g.demultiplexed[sample]['index2'] = os.path.join(g.analysis_folder, 'demultiplexed', sample + '.i2.fastq')
                    if not os.path.isfile(g.demultiplexed[sample]['read1']):
                        print ("Can't find ",g.demultiplexed[sample]['read1'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['read2']):
                        print ("Can't find ",g.demultiplexed[sample]['read2'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['index1']):
                        print ("Can't find ",g.demultiplexed[sample]['index1'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['index2']):
                        print ("Can't find ",g.demultiplexed[sample]['index2'])
                        exit()

                    sample = "control_"+sample    
                    g.demultiplexed[sample] = {}
                    g.demultiplexed[sample]['read1'] = os.path.join(g.analysis_folder, 'demultiplexed', sample + '.r1.fastq')
                    g.demultiplexed[sample]['read2'] = os.path.join(g.analysis_folder, 'demultiplexed', sample + '.r2.fastq')
                    g.demultiplexed[sample]['index1'] = os.path.join(g.analysis_folder, 'demultiplexed', sample + '.i1.fastq')
                    g.demultiplexed[sample]['index2'] = os.path.join(g.analysis_folder, 'demultiplexed', sample + '.i2.fastq')
                    if not os.path.isfile(g.demultiplexed[sample]['read1']):
                        print ("Can't find ",g.demultiplexed[sample]['read1'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['read2']):
                        print ("Can't find ",g.demultiplexed[sample]['read2'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['index1']):
                        print ("Can't find ",g.demultiplexed[sample]['index1'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['index2']):
                        print ("Can't find ",g.demultiplexed[sample]['index2'])
                        exit()


                g.umitag()
                g.consolidate()
                g.alignReads()
                g.identifyOfftargetSites()
                g.visualize()

            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        else:
            g = GuideSeq()
            g.parseManifest(args.manifest)
            g.demultiplex()
            g.umitag()
            g.consolidate()
            g.alignReads()
            g.identifyOfftargetSites()
            g.visualize()

    elif args.command == 'demultiplex':
        """
        Run just the demultiplex step given the manifest
        """
        g = GuideSeq()
        g.parseManifest(args.manifest)
        g.demultiplex()
        g.umitag()
        g.consolidate()
 

    elif args.command == 'align':
        """
        Run just the alignment step
        python guideseq/guideseq.py align --bwa bwa --read1 test/data/consolidated/EMX1.r1.consolidated.fastq --read2 test/data/consolidated/EMX1.r2.consolidated.fastq --genome /Volumes/Media/hg38/hg38.fa --outfolder test/output/
        """
        g = GuideSeq()
        g.parseManifest(args.manifest)
        g.consolidated = {sample: {}}
        g.consolidated[sample]['read1'] = args.read1
        g.consolidated[sample]['read2'] = args.read2
        g.alignReads()

    elif args.command == 'identify':
        """
        Run just the identify step
        python guideseq/guideseq.py identify --genome /Volumes/Media/hg38/hg38.fa --aligned test/output/aligned/EMX1.sam --outfolder test/output/ --target_sequence GAGTCCGAGCAGAAGAAGAANGG
        """
        if 'description' in args:
            description = args.description
        else:
            description = ''

        if 'max_score' in args:
            max_score = args.max_score
        else:
            max_score = 7

        if 'window_size' in args:
            window_size = args.window_size
        else:
            window_size = 25

        g = GuideSeq()
        g.parseManifest(args.manifest)
        sample = os.path.basename(args.aligned).split('.')[0]
        g.samples = {sample: {'description': description, 'target': args.target_sequence}}
        g.aligned = {sample: args.aligned}
        g.identifyOfftargetSites()


    elif args.command == 'visualize':
        """
        Run just the visualize step
        """
        g = GuideSeq()
        g.parseManifest(args.manifest)
        g.identified = {}
        g.identified[sample] = args.infile
        g.visualize()


if __name__ == '__main__':
    main()
