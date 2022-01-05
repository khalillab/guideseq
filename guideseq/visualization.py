from __future__ import print_function
import svgwrite
import sys
import os
import logging

logger = logging.getLogger('root')
logger.propagate = False

boxWidth = 10
box_size = 15
v_spacing = 3

# colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', 'R': '#B3B3B3', '-': '#FFFFFF'}
colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', 'R': '#B3B3B3', '-': '#B3B3B3'}
for c in ['Y','S','W','K','M','B','D','H','V','.']:
    colors[c] = "#B3B3B3"

def parseSitesFile(infile):
    sample_names = []
    offtargets = []
    total_seq = 0
    target_seqs = []
    if os.path.isdir(infile):
        files = os.listdir(infile)
    else:
        files = [infile]
    for input_file in files:
        if not os.path.isfile(os.path.join(infile, input_file)) or 'DS_Store' in input_file:
            continue
        with open(os.path.join(infile, input_file), 'r') as f:
            f.readline()
            for line in f:
                line = line.rstrip('\n')
                line_items = line.split('\t')
                offtarget_reads = line_items[11]
                no_bulge_offtarget_sequence = line_items[24]
                bulge_offtarget_sequence = line_items[29]
                target_seq = line_items[40]
                if target_seq not in target_seqs:
                    target_seqs.append(target_seq)
                realigned_target_seq = line_items[41]
                sample_name = line_items[42]
                pam = line_items[43]
                if sample_name not in sample_names:
                    sample_names.append(sample_name)

                if no_bulge_offtarget_sequence != '' or bulge_offtarget_sequence != '':
                    if no_bulge_offtarget_sequence:
                        total_seq += 1
                    if bulge_offtarget_sequence:
                        total_seq += 1
                    offtargets.append({'seq': no_bulge_offtarget_sequence.strip(),
                                       'bulged_seq': bulge_offtarget_sequence.strip(),
                                       'reads': [int(offtarget_reads.strip())],
                                       'target_seq': target_seq.strip(),
                                       'realigned_target_seq': realigned_target_seq.strip(),
                                       'sample_name': [sample_name.strip()],
                                       'pam': pam.strip(),
                                       'other': []
                                       })
    to_remove = []
    offtargets = sorted(offtargets, key=lambda x: x['reads'][0], reverse=True)
    sample_names = sorted(sample_names)
    seqs = {}
    
    for i, offtarget in enumerate(offtargets):
        if offtarget['seq'] in seqs.keys():
            if not offtarget['sample_name'][0] in offtargets[seqs[offtarget['seq']]]['sample_name']:
                offtargets[seqs[offtarget['seq']]]['reads'].append(offtarget['reads'][0])
                offtargets[seqs[offtarget['seq']]]['sample_name'].append(offtarget['sample_name'][0])
                offtargets[seqs[offtarget['seq']]]['other'].append(offtarget)
            to_remove.append(offtarget)
        else:
            seqs[offtarget['seq']] = i
        """
        elif offtarget['pam'] == 'NNNNNN':
            if offtarget['seq'][3:23] in seqs.keys():
                if not offtarget['sample_name'][0] in offtargets[seqs[offtarget['seq'][3:23]]]['sample_name']:
                    offtargets[seqs[offtarget['seq'][3:23]]]['reads'].append(offtarget['reads'][0])
                    offtargets[seqs[offtarget['seq'][3:23]]]['sample_name'].append(offtarget['sample_name'][0])
                    offtargets[seqs[offtarget['seq'][3:23]]]['other'].append(offtarget)
                to_remove.append(offtarget)
            else:
                seqs[offtarget['seq']] = i
                seqs[offtarget['seq'][3:23]] = i
        elif offtarget['seq'][0:20] in seqs.keys():
            if not offtarget['sample_name'][0] in offtargets[seqs[offtarget['seq'][0:20]]]['sample_name']:
                offtargets[seqs[offtarget['seq'][0:20]]]['reads'].append(offtarget['reads'][0])
                offtargets[seqs[offtarget['seq'][0:20]]]['sample_name'].append(offtarget['sample_name'][0])
                offtargets[seqs[offtarget['seq'][0:20]]]['other'].append(offtarget)
            to_remove.append(offtarget)
            """
    for offtarget in to_remove:
        try:
            offtargets.remove(offtarget)
        except ValueError:
            continue
    print(offtargets)
    return offtargets, list(sorted(target_seqs, key = lambda x: len(x))), total_seq, sample_names

# 3/6/2020
def check_mismatch(a,b):
    from Bio.Data import IUPACData
    try:
        dna_dict = IUPACData.ambiguous_dna_values
        set_a = dna_dict[a.upper()]
        set_b = dna_dict[b.upper()]
        overlap = list(set(list(set_a)).intersection(list(set_b)))
    except:
        return True
    if len(overlap) == 0:
        return True
    else:
        return False

def draw_aligned_sequence_row(line_number, j, seq_offset, target_seq_to_use, realigned_target_seq, no_bulge_offtarget_sequence, bulge_offtarget_sequence, pam, target_seqs, dwg, second=False):
    x_offset = 20
    y_offset = 50

    if second:
        print(target_seq_to_use)

    if no_bulge_offtarget_sequence != '':
        k = 0
        line_number += 1
        y = y_offset + line_number * box_size
        align_offset = 0
        for i, (c, r) in enumerate(zip(no_bulge_offtarget_sequence, target_seq_to_use)):
            x = x_offset + k * box_size + box_size * align_offset + seq_offset * 110 * 5 + seq_offset * (len(target_seqs[0]) * box_size + 16)
            if r == '-':
                if 0 < k < len(target_seq_to_use):
                    x = x_offset + (k - 0.25) * box_size + align_offset + seq_offset * 110 * 5 + seq_offset * (len(target_seqs[0]) * box_size + 16)
                    # if i < PAM_index:
                    if check_mismatch(c,r):
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
                    else:
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill="#FFFFFF"))
                    dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
            elif c == r:
                dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
                k += 1
            elif r == 'N':
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                k += 1
            else:
                if check_mismatch(c,r):
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                else:
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill="#FFFFFF"))
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                k += 1
    if bulge_offtarget_sequence != '':
        k = 0
        line_number += 1
        y = y_offset + line_number * box_size
        align_offset = 0
        for i, (c, r) in enumerate(zip(bulge_offtarget_sequence, realigned_target_seq)):
            x = x_offset + k * box_size + box_size * align_offset + seq_offset * 110 * 5 + seq_offset * (len(target_seqs[0]) * box_size + 16)
            if r == '-':
                if 0 < k < len(realigned_target_seq):
                    x = x_offset + (k - 0.25) * box_size + align_offset + seq_offset * 110 * 5 + seq_offset * (len(target_seqs[0]) * box_size + 16)
                    if check_mismatch(c,r):
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
                    else:
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill="#FFFFFF"))
                    dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
            elif c == r:
                dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
                k += 1
            elif r == 'N':
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                k += 1
            else:
                if check_mismatch(c,r):
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                else:
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill="#FFFFFF"))
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                k += 1
    return line_number

def visualizeOfftargets(infile, outfile, title, PAMS):

    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get offtargets array from file
    offtargets, target_seqs, total_seq, sample_names = parseSitesFile(infile)

    # Initiate canvas
    dwg = svgwrite.Drawing(outfile + '.svg', profile='full', size=(u'100%', 100 + total_seq*(box_size + 1)))

    if title is not None:
        # Define top and left margins
        x_offset = 20
        y_offset = 50
        dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
    else:
        # Define top and left margins
        x_offset = 20
        y_offset = 20

    # Draw ticks
    # tick_locations = [1, len(target_seq)]  # limits
    # if target_seq.index('N') > len(target_seq)/2:  # PAM on the right end
        # tick_locations += range(len(target_seq) + 1)[::10][1:]  # intermediate values
        # tick_locations += range(len(target_seq) + 1)[len(target_seq) - 2: len(target_seq)]  # complementing PAM
        # tick_locations.sort()
        # tick_legend = [str(x) for x in tick_locations[:-3][::-1]] + ['P', 'A', 'M']
    # else:
        # tick_locations += [range(3, len(target_seq) + 1)[::10][1]]
        # tick_locations += range(2, 5)
        # tick_locations.sort()
        # tick_legend = ['P', 'A', 'M'] + [str(x) for x in [str(x-3) for x in tick_locations[3:]]]
    ## Assume PAM is on the right end
    print(target_seqs)
    for ts_iterator, target_seq in enumerate(target_seqs):
        tick_locations = []
        tick_legend = []
        try:
            PAM_index = target_seq.index(PAMS[0])
        except ValueError:
            PAM_index = target_seq.index(PAMS[1])
        count = 0
        for i in range(PAM_index,0,-1):
            count = count+1
            if count % 10 == 0:
                tick_legend.append(count)
                tick_locations.append(i)
        tick_legend+=['P', 'A', 'M']+['-']*(len(PAMS[ts_iterator])-3)
        tick_locations+=range(PAM_index+1,len(target_seq)+1)
        


        for x,y in zip(tick_locations, tick_legend):
            dwg.add(dwg.text(y, insert=(x_offset + (x - 1) * box_size + 2 + ts_iterator * 110 * 5 + ts_iterator * (len(target_seqs[0]) * box_size + 16), y_offset - 2), style="font-size:10px; font-family:Courier"))

        # Draw reference sequence row
        for i, c in enumerate(target_seq):
            y = y_offset
            x = x_offset + i * box_size
            if i < PAM_index:
                dwg.add(dwg.rect((x + ts_iterator * 110 * 5 + ts_iterator * (len(target_seqs[0]) * box_size + 16), y), (box_size, box_size), fill=colors[c]))
            else:
                dwg.add(dwg.rect((x + ts_iterator * 110 * 5 + ts_iterator * (len(target_seqs[0]) * box_size + 16), y), (box_size, box_size), fill="#B3B3B3"))
            dwg.add(dwg.text(c, insert=(x + 3 + ts_iterator * 110 * 5 + ts_iterator * (len(target_seqs[0]) * box_size + 16), y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))

    if sample_names is not None:
        scount = 0
        for sample_name in sample_names:
            dwg.add(dwg.text(sample_name, insert=(x_offset + box_size * len(target_seqs[0]) + 16 + scount * 110, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
            scount += 1
    else:
        dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(target_seqs[0]) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))

    # Draw aligned sequence rows
    y_offset += 1  # leave some extra space after the reference row
    line_number = 0  # keep track of plotted sequences
    last_line_number = 0
    for j, seq in enumerate(offtargets):
        seq_offset = 0
        if offtargets[j]['target_seq'] == target_seqs[0]:
            target_seq_to_use = target_seqs[0]
        else:
            target_seq_to_use = target_seqs[1]
            seq_offset = 1
        realigned_target_seq = offtargets[j]['realigned_target_seq']
        no_bulge_offtarget_sequence = offtargets[j]['seq']
        bulge_offtarget_sequence = offtargets[j]['bulged_seq']
        pam = offtargets[j]['pam']

        sub = 1
        line_number = draw_aligned_sequence_row(line_number, j, seq_offset, target_seq_to_use, realigned_target_seq, no_bulge_offtarget_sequence, bulge_offtarget_sequence, pam, target_seqs, dwg)
        if line_number - last_line_number > 1:
            sub = 2

        last_line_number = line_number

        for other in offtargets[j]['other']:
            seq_offset = 0
            if other['target_seq'] == target_seqs[0]:
                target_seq_to_use = target_seqs[0]
            else:
                target_seq_to_use = target_seqs[1]
                seq_offset = 1
            realigned_target_seq = other['realigned_target_seq']
            no_bulge_offtarget_sequence = other['seq']
            #bulge_offtarget_sequence = other['bulged_seq']
            bulge_offtarget_sequence = ''
            pam = other['pam']

            draw_aligned_sequence_row(line_number-sub, j, seq_offset, target_seq_to_use, realigned_target_seq, no_bulge_offtarget_sequence, bulge_offtarget_sequence, pam, target_seqs, dwg)


        if no_bulge_offtarget_sequence == '' or bulge_offtarget_sequence == '':
            for k,sample_name in enumerate(offtargets[j]['sample_name']):
                for s in sample_names:
                    if s in offtargets[j]['sample_name']:
                        reads_text = dwg.text(str(seq['reads'][k]), insert=(box_size * (len(target_seqs[0]) + 1) + 20 + sample_names.index(sample_name)*110, y_offset + box_size * (line_number + 2) - 2), fill='black', style="font-size:15px; font-family:Courier")
                        dwg.add(reads_text)
                    else:
                        reads_text = dwg.text('--', insert=(box_size * (len(target_seqs[0]) + 1) + 20 + sample_names.index(s)*110, y_offset + box_size * (line_number + 2) - 2), fill='black', style="font-size:15px; font-family:Courier")
                        dwg.add(reads_text)
        else:
            for k, sample_name in enumerate(offtargets[j]['sample_name']):
                for s in sample_names:
                    if s in offtargets[j]['sample_name']:
                        reads_text = dwg.text(str(seq['reads'][k]), insert=(box_size * (len(target_seqs[0]) + 1) + 20 + sample_names.index(sample_name)*110, y_offset + box_size * (line_number + 1) + 5), fill='black', style="font-size:15px; font-family:Courier")
                        dwg.add(reads_text)
                    else:
                        reads_text = dwg.text('--', insert=(box_size * (len(target_seqs[0]) + 1) + 20 + sample_names.index(s)*110, y_offset + box_size * (line_number + 1) + 5), fill='black', style="font-size:15px; font-family:Courier")
                        dwg.add(reads_text)

            reads_text02 = dwg.text(u"\u007D", insert=(box_size * (len(target_seqs[seq_offset]) + 1) + 7 + seq_offset * (len(target_seqs[0]) * box_size + 16) + seq_offset * 110 * 5, y_offset + box_size * (line_number + 1) + 5), fill='black', style="font-size:23px; font-family:Courier")
            dwg.add(reads_text02)
    dwg.save()


def main():
    if len(sys.argv) >= 3:
        if len(sys.argv) == 5:
            title = sys.argv[3]
            PAM = sys.argv[4]
        else:
            title = None
            PAM="NGG"
        visualizeOfftargets(sys.argv[1], sys.argv[2], title=title,PAM=PAM)
    else:
        print('Usage: python visualization.py INFILE OUTFILE [TITLE] [PAM]')

if __name__ == '__main__':
    main()
