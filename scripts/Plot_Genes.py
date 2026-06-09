#!/bin/python

import math
import cairo
import argparse
import sys
import gzip
import copy

def get_arguments():
    """get the arguments"""

    parser = argparse.ArgumentParser(
        description="Plot gene of interest with cmh values after FDR")
    parser.add_argument(
        "-g", "--gtf", help="Input GTF file", required=True)
    parser.add_argument(
        "-c", "--cmh", help="cmh file", required=True)
    parser.add_argument(
        "-G", "--gap_file", help="gaps file", default=None)
    parser.add_argument(
        "-l", "--leftbound", help="lefbound", default=0, type=int)
    parser.add_argument(
        "-r", "--rightbound", help="rightbound", default=10000000,type=int)
    parser.add_argument(
        "-C", "--chrom", help="name of chromosome to plot", default="HiC_scaffold_12")

    args = parser.parse_args()
    gtf = args.gtf
    cmh = args.cmh
    gap_file = args.gap_file
    leftbound = args.leftbound
    rightbound = args.rightbound
    chrom = args.chrom

    assert leftbound < rightbound, "Leftbound must be smaller than rightbound"

    return gtf, cmh, gap_file, leftbound, rightbound, chrom

class SNP:
    def __init__(self, pos, value):
        self.pos = pos
        self.value = value

class FEATURE:
    def __init__(self, left_pos, right_pos, direction, color):
        self.lpos = left_pos
        self.rpos = right_pos
        self.direction = direction
        self.color = color

#################
### functions ###
#################

######### parsing functions ############
def parse_cmh(cmh_file, leftbound, rightbound):

    snps = list()
    max_cmh = 0

    if str(cmh_file).endswith(".gz"):
        fh = gzip.open(cmh_file, "rt")
    else:
        fh = open(cmh_file, 'r')

    for line in fh:
        fields = line.strip('\n').split('\t')
        pos = int(fields[1])
        if leftbound <= pos <= rightbound:
            value = float(fields[9])
            if value > max_cmh:
                max_cmh = value
            snps.append(SNP(pos, value))
    fh.close()

    return snps, max_cmh

def get_gene_name(gene_info):
    while (gene_info.startswith("gene_id") == False):
        gene_info = gene_info[1:]
    gene_name = gene_info.split(';')[0].split(' ')[1].replace('"', '')
    return gene_name

def get_color(num):

    colors = {"cinnabar": [0.89, 0.26, 0.20],
              "maroon": [0.50, 0, 0],
              "crimson": [0.86, 0.08, 0.24],      
              "indigo": [0.29, 0, 0.51],
              "lavender": [0.71, 0.49, 0.86],
              "mango": [0.99, 0.75, 0.01], 
              "peach": [1, 0.90, 0.71],
              "periwinkle": [0.80, 0.80, 1],
              "sapphire": [0.06, 0.32, 0.73],
              "cyan": [0,1,1],    
            }

    colors_list = list(colors.keys())

    color = colors_list[num]

    return colors[color]


def parse_gtf(gtf_file, leftbound, rightbound, chrom):

    if str(gtf_file).endswith(".gz"):
        fh = gzip.open(gtf_file, "rt")
    else:
        fh = open(gtf_file, 'r')

    exon_count = {}
    genes_found = list()
    gene_bounds = {}

    found_gene = False

    exons = {}

    for line in fh:
        fields = line.strip('\n').split('\t')
        if fields[0] == chrom:
            if fields[2] == "gene":
                left_pos = int(fields[3])
                right_pos = int(fields[4])
                if leftbound <= left_pos <= rightbound:
                    if right_pos <= rightbound:
                        pass
                    else:
                        right_pos = rightbound
                    found_gene = True
                    gene_name = get_gene_name(fields[8])
                    exon_count[gene_name] = 0
                    gene_bounds[gene_name] = [left_pos, right_pos]
                elif leftbound <= right_pos <= rightbound:
                    left_pos = leftbound
                    found_gene = True
                    gene_name = get_gene_name(fields[8])
                    exon_count[gene_name] = 0
                    gene_bounds[gene_name] = [left_pos, right_pos]
                else:
                    found_gene = False
            elif found_gene:
                if fields[2] == "exon":
                    left_pos = int(fields[3])
                    right_pos = int(fields[4])
                    gene_name = get_gene_name(fields[8])
                    gene_bd = gene_bounds[gene_name]
                    b1 = gene_bd[0]
                    b2 = gene_bd[1]
                    if right_pos < b1 or left_pos > b2:
                        continue
                    elif left_pos < b1 and right_pos < b2:
                        left_pos = b1
                    elif left_pos < b2 and right_pos > b2:
                        right_pos = b2
                    left_pos = left_pos - leftbound
                    right_pos = right_pos - leftbound
                    if left_pos < 0 or right_pos < 0:
                        continue
                    if gene_name not in genes_found:
                        genes_found.append(gene_name)
                    num_genes = len(genes_found)
                    num = num_genes % 5
                    gene_color = get_color(num)
                    direction = fields[6]
                    exon_count[gene_name] += 1
                    ex = FEATURE(left_pos, right_pos, direction, gene_color)
                    cnt = exon_count[gene_name]
                    exons[f"{gene_name}_E{cnt}"] = ex
    # print gene names to console
    genes = ' '.join(genes_found)
    print("Found the following genes: " + genes)

    return exons

def parse_gap_file(chrom, gap_file, leftbound, rightbound):

    gaps = list()
    gap_color = [48/255, 48/255, 48/255]

    if str(gap_file).endswith(".gz"):
        fh = gzip.open(gap_file, "rt")
    else:
        fh = open(gap_file, 'r')

    for line in fh:
        fields = line.strip('\n').split('\t')
        if fields[0] == chrom:
            in_region = False
            left_pos = int(fields[1])
            right_pos = int(fields[2])
            if leftbound <= left_pos <= rightbound:
                in_region = True
                if right_pos <= rightbound:
                    pass
                else:
                    right_pos = rightbound
            elif leftbound <= right_pos <= rightbound:
                left_pos = leftbound
                in_region = True
            if in_region:
                left_pos = left_pos - leftbound
                right_pos = right_pos - leftbound
                gap = FEATURE(left_pos, right_pos, '+', gap_color)
                gaps.append(gap)

    return gaps            

########## helper functions #############

def get_y_spacing(pdf_height, y_offset, max_cmh):
    """adjust pdf height for # of seqs"""

    max_cmh = math.ceil(max_cmh)
    ceiling = copy.copy(max_cmh)

    while (ceiling % 5 != 0):
        ceiling += 1

    num_ticks = (ceiling / 5) + 1 # include tick for 0

    space = pdf_height - y_offset

    interval = space / num_ticks

    return interval, ceiling

########## plotting functions #############


def plot_chrom(cairo_context, seq_length, seq_id, allocated_pos, seq_height, x_pos, leftbound, rightbound):
    """ plot sequence to scale to the longest sequence"""

    sequence_id = seq_id

    # starting positions
    x1 = x_pos  # offset from left by x position
    y1 = seq_height  # y / vertical offset

    # determine new x
    x2 = x_pos + (allocated_pos * seq_length)  # scale to longest seq
    y2 = y1

    cairo_context.move_to(x1, y1)
    cairo_context.line_to(x2, y2)

    # create seq color (i.e., light grey)
    seq_color = [0.7, 0.7, 0.7, 1]
    cairo_context.set_source_rgba(
        seq_color[0], seq_color[1], seq_color[2], seq_color[3]
    )
    cairo_context.set_line_width(50)
    cairo_context.stroke()

    # first horizontal line
    cairo_context.set_source_rgb(0, 0, 0)  # black
    cairo_context.move_to(x1, y1 + 23)
    cairo_context.line_to(x2, y1 + 23)
    cairo_context.set_line_width(5)
    cairo_context.stroke()

    # second horizontal line
    cairo_context.set_source_rgb(0, 0, 0)
    cairo_context.move_to(x1, y1 - 23)
    cairo_context.line_to(x2, y1 - 23)
    cairo_context.set_line_width(5)
    cairo_context.stroke()

    # first veritcal line
    cairo_context.set_source_rgb(0, 0, 0)
    cairo_context.move_to(x1, y1 + 23)
    cairo_context.line_to(x1, y1 - 23)
    cairo_context.set_line_width(5)
    cairo_context.stroke()

    # second vertical line
    cairo_context.set_source_rgb(0, 0, 0)
    cairo_context.move_to(x2, y1 + 23)
    cairo_context.line_to(x2, y1 - 23)
    cairo_context.set_line_width(5)
    cairo_context.stroke()

    # add sequence name based if align_pos to top
    x, y = 100, seq_height + 160

    cairo_context.set_source_rgb(0, 0, 0)
    cairo_context.set_font_size(75)
    cairo_context.select_font_face(
        "Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL
    )
    cairo_context.move_to(x, y)
    left = "{:,}".format(leftbound)
    right = "{:,}".format(rightbound)
    seq_id = f"{sequence_id}: {left}-{right}"
    cairo_context.show_text(seq_id)
    cairo_context.stroke()


def create_ticks(seq_length):
    """Find the appropriate marker lengths for tick marks"""

    marker_sets = {
        1e3: ["bp", 100],
        5e3: ["K", 1000],
        7.5e3: ["K", 1500],
        1e4: ["K", 2000],
        2.5e4: ["K", 5000],
        5e4: ["K", 10000],
        9e4: ["K", 15000],
        1e5: ["K", 20000],
        5e5: ["K", 25000],
        9e5: ["K", 30000],
        1e6: ["M", 50000],
        5e6: ["M", 70000],
        9e6: ["M", 100000],
        1e7: ["M", 150000],
        3e7: ["M", 200000],
        5e7: ["M", 500000],
        7e7: ["M", 700000],
        9e7: ["M", 900000],
        1e8: ["M", 1000000],
        5e8: ["M", 1500000],
        9e8: ["M", 3000000],
        1e9: ["M", 5000000],
              }

    multiplier = None

    for marker_length in marker_sets.keys():
        if seq_length <= marker_length:
            fields = marker_sets[marker_length]
            multiplier = fields[1]
            break

    # multiplier will be used to show increments in bp length
    if multiplier != None:
        tick_ct = math.floor(seq_length/multiplier)
    else:
        sys.exit(f"Software was not prepared for sequence/segment of length {seq_length}")

    return tick_ct, multiplier


def draw_ticks(cairo_context, y_coordinate, seq_length, image_width, section, x_pos, leftbound):
    """add tick marks to the reference sequence"""

    # need to offset from the chrom lines
    mark_ypos1 = y_coordinate + (section * 0.4)
    mark_ypos2 = y_coordinate + (section * 1.08)
    label_ypos = y_coordinate + (section * 1.55)

    tick_ct, multiplier = create_ticks(seq_length)

    # adjust for just ticks
    allocated_pos = ((image_width / seq_length) * \
        seq_length) / (seq_length/multiplier)


    for i in range(1, tick_ct + 1):
        marker_num = i * multiplier  # interval scheme
        marker_num = marker_num + leftbound #adjust to ref
        # for rounding
        if marker_num < 1e3:
            marker_num = marker_num/100
            marker = "bp"
        elif 1e3 < marker_num < 1e6:
            marker_num = marker_num/1e3
            marker = "K"
        elif marker_num > 1e6:
            marker_num = marker_num/1e6
            marker = "M" 
        marker_num = round(marker_num, 2)
        marker_label = f"{str(marker_num)}{marker}"
        marker_pos = x_pos + (allocated_pos * i)
        cairo_context.move_to(marker_pos, mark_ypos1)
        cairo_context.line_to(marker_pos, mark_ypos2)
        cairo_context.set_source_rgba(0, 0, 0, 1)  # black
        cairo_context.set_line_width(section * 0.05)
        cairo_context.stroke()
        # marker labels
        cairo_context.set_source_rgba(0, 0, 0, 1)
        cairo_context.set_font_size(45)
        cairo_context.move_to(marker_pos - (section * 0.45), label_ypos)
        cairo_context.show_text(marker_label)
        cairo_context.stroke()


def plot_label(cairo_context, y, theta, label, midpoint):
    """plot clump labels"""

    new_x = midpoint
    new_y = y - 90 # 63
    cairo_context.save()
    # create labels
    cairo_context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL,
                                   cairo.FONT_WEIGHT_NORMAL)
    cairo_context.set_source_rgb(0, 0, 0)  # black
    cairo_context.set_font_size(25)
    cairo_context.move_to(new_x, new_y)
    cairo_context.rotate(theta)
    cairo_context.show_text(label)
    cairo_context.stroke()
    cairo_context.restore()

def draw_y_axis(cairo_context, x_pos, pdf_height, y_offset, spacing, ceiling, max_cmh):

    num_ticks = int(ceiling / 5) + 1

    start_y = pdf_height - (y_offset * 1.75)

    top_y = start_y - ((num_ticks - 1) * spacing)

    # draw vertical line
    cairo_context.move_to(x_pos, start_y)
    cairo_context.line_to(x_pos, top_y)
    cairo_context.set_source_rgb(0, 0, 0) 
    cairo_context.set_line_width(8)
    cairo_context.stroke()

    # now for the ticks
    for n in range(num_ticks):
        y_pos = start_y - (n * spacing)
        cairo_context.move_to(x_pos, y_pos)
        cairo_context.line_to(x_pos + 50, y_pos)
        cairo_context.set_source_rgb(0, 0, 0) 
        cairo_context.set_line_width(8)
        cairo_context.stroke()
        # label
        cairo_context.set_source_rgb(0, 0, 0)
        cairo_context.set_font_size(75)
        cairo_context.select_font_face(
            "Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL
        )
        cairo_context.move_to(x_pos - 100, y_pos)
        cairo_context.show_text(str(n * 5))
        cairo_context.stroke()

    # now add a y-axis label
    y_mid = (start_y - top_y) * 0.65
    ylab = "-log(qvalue)"
    theta = math.radians(270)
    cairo_context.set_source_rgb(0, 0, 0)
    cairo_context.set_font_size(90)
    cairo_context.select_font_face(
        "Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL
    )
    cairo_context.move_to(x_pos - 200, y_mid)
    cairo_context.rotate(theta)
    cairo_context.show_text(ylab)
    cairo_context.stroke()


def plot_exons(cairo_context, y_coordinate, allocated_pos, exons, x_pos):
    """plot exons"""

    for ex, exon in exons.items():
        fill_in_matches(cairo_context, y_coordinate, x_pos,
                        allocated_pos, exon.lpos, exon.rpos, exon.color)
        # get direction to see if clump is on forward or reverse
        direction = exon.direction
        if direction == "+":
            yLabel = y_coordinate + 60
            angle = -45
        elif direction == "-":
            yLabel = y_coordinate + 105
            angle = 45
        # add label
        label = ex
        # convert to radians
        theta = math.radians(angle)
        # get the mid point of the text if i were not rotated
        start_pos = x_pos + (exon.lpos * allocated_pos)
        end_pos = x_pos + (exon.rpos * allocated_pos)
        midpoint = (end_pos + start_pos) / 2
        # create labels
        plot_label(cairo_context, yLabel, theta, label, midpoint)


def fill_in_matches(cairo_context, y, x_pos, allocated_pos, lpos, rpos, color):
    """fill in seq with clump positions"""

    for i in range(lpos, rpos + 1):
        site = x_pos + (i * allocated_pos)  # scale
        cairo_context.move_to(site, y - 20)
        cairo_context.line_to(site, y + 20)
        cairo_context.set_source_rgb(color[0], color[1], color[2]) 
        cairo_context.set_line_width(5)
        cairo_context.stroke()

def plot_snps(cairo_context, snps, y_pos, x_pos, allocated_pos, leftbound, color):

    for snp in snps:
        pos = x_pos + ((snp.pos - leftbound) * allocated_pos)
        cairo_context.move_to(pos, y_pos - 20)
        cairo_context.line_to(pos, y_pos + 20)
        cairo_context.set_source_rgb(color[0], color[1], color[2]) 
        cairo_context.set_line_width(5)
        cairo_context.stroke()

def plot_cmh(cairo_context, snps, pdf_height, y_offset, x_start, allocated_pos, leftbound, ceiling, spacing):

    num_ticks = int(ceiling / 5) + 1

    allocated_pos_y = pdf_height - (y_offset * 1.75)
    # allocated_pos_y = allocated_pos_y - (((num_ticks - 1) * spacing))

    for snp in snps:
        x_pos = x_start + ((snp.pos - leftbound) * allocated_pos)
        y_pos = (1 - (snp.value / ceiling)) * allocated_pos_y
        cairo_context.arc(x_pos, y_pos, 15, 0, 2*math.pi)
        cairo_context.close_path()
        cairo_context.set_source_rgb(0, 0, 0)
        cairo_context.fill()


def draw_image(exons, leftbound, rightbound, snps, gaps, chrom, max_cmh):
    """plot clumps using plotting functions"""

    pdf_width = 5000  # 7500
    pdf_height = 3000
    y_offset = 250
    chrom_y_pos = pdf_height - y_offset
    chrom_x_pos = 450
    y_axix_x_pos = 400
    section = pdf_height * 0.02
    
    spacing, ceiling = get_y_spacing(pdf_height, y_offset, max_cmh)
    seq_len = rightbound - leftbound

    ims = cairo.PDFSurface(
        f"{chrom}-{leftbound}-{rightbound}.pdf", pdf_width, pdf_height)
    cairo_context = cairo.Context(ims)

    image_width = pdf_width - (chrom_x_pos * 1.5)  # offset by starting positions on both sides

    # create allocated posititions for sites
    allocated_pos = image_width / seq_len

    # plot seqs
    plot_chrom(cairo_context, seq_len, chrom, allocated_pos, chrom_y_pos, chrom_x_pos, leftbound, rightbound)
    draw_ticks(cairo_context, chrom_y_pos, seq_len, image_width, section, chrom_x_pos, leftbound)
    plot_exons(cairo_context, chrom_y_pos, allocated_pos, exons, chrom_x_pos)
    if gaps != None and len(gaps) > 0:
        for gap in gaps:
            fill_in_matches(cairo_context, chrom_y_pos, chrom_x_pos, allocated_pos, gap.lpos, gap.rpos, gap.color)
    # plot_snps(cairo_context, snps, chrom_y_pos, chrom_x_pos, allocated_pos, leftbound, snp_colors)
    plot_cmh(cairo_context, snps, pdf_height, y_offset, chrom_x_pos, allocated_pos, leftbound, ceiling, spacing)
    draw_y_axis(cairo_context, y_axix_x_pos, pdf_height, y_offset, spacing, ceiling, max_cmh)
    
    
    # finish plot
    cairo_context.show_page()


def main():
    """Parse, plot, and pray for good results"""

    # get arguments
    gtf, cmh, gap_file, leftbound, rightbound, chrom = get_arguments()

    # parse inputs
    exons = parse_gtf(gtf, leftbound, rightbound, chrom)

    snps, max_cmh = parse_cmh(cmh, leftbound, rightbound)

    if gap_file != None:
        gaps = parse_gap_file(chrom, gap_file, leftbound, rightbound)
    else:
        gaps = None

    # run the entire script
    draw_image(exons, leftbound, rightbound, snps, gaps, chrom, max_cmh)

if __name__ == "__main__":
    main()
