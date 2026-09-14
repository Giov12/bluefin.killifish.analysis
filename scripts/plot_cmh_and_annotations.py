#!/bin/python

import math
import cairo
import argparse
import sys
import gzip
import os
import copy

# create some global variables
cmh        = ''
gtf        = ''
chrom      = ''
gap_file   = ''
leftbound  = 0
rightbound = 10_000_000

def get_arguments() -> int:
    """get the arguments"""

    global cmh, gtf, chrom, gap_file, leftbound, rightbound

    parser = argparse.ArgumentParser(description="Plot region of interest with cmh values after FDR")
    parser.add_argument("-g", "--gtf", help="Gene annotations in GTF/GFF3 format", required=True)
    parser.add_argument("-c", "--cmh", help="cmh file", required=True)
    parser.add_argument("-G", "--gap_file", help="gaps file generated via Klumpy[optional]", default='')
    parser.add_argument("-l", "--leftbound", help="lefbound", default=leftbound, type=int)
    parser.add_argument("-r", "--rightbound", help="rightbound", default=rightbound,type=int)
    parser.add_argument("-C", "--chrom", help="name of chromosome to plot", default=chrom, required=True)

    args       = parser.parse_args()
    gtf        = args.gtf
    cmh        = args.cmh
    gap_file   = args.gap_file
    leftbound  = args.leftbound
    rightbound = args.rightbound
    chrom      = args.chrom

    assert os.path.isfile(cmh), f"Could not locate {cmh}"
    assert os.path.isfile(gtf), f"Could not locate {gtf}"
    assert chrom != '', "Chromosome name needs to be specified"
    assert leftbound > 0, "leftbound must be greater than 0"
    assert leftbound < rightbound, "Leftbound must be smaller than rightbound"

    if (gap_file != ''):
        assert os.path.isfile(gap_file), f"Could not locate {gap_file}"

    return 0

class SNP:
    def __init__(self, pos: int, value: float):
        self.pos = pos
        self.value = value

class FEATURE:
    def __init__(self, left_pos: int, right_pos: int, direction: str, color: list[float]):
        self.lpos      = left_pos
        self.rpos      = right_pos
        self.direction = direction
        self.color     = color

#################
### functions ###
#################

######### parsing functions ############
def parse_cmh() -> tuple[list[SNP], float]:
    """get the SNPs"""

    global cmh, leftbound, rightbound

    snps    = list()
    max_cmh = 0

    if (cmh.endswith(".gz")):
        fh = gzip.open(cmh, "rt")
    else:
        fh = open(cmh, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.strip('\n').split('\t')
        pos    = int(fields[1])
        if (leftbound <= pos <= rightbound):
            value = float(fields[2])
            if (value > max_cmh):
                max_cmh = value
            snps.append(SNP(pos, value))
    fh.close()

    # max was found to be 31.7304182677135

    return snps, max_cmh

def get_gene_id(attrb: str) -> str:
    """return the gene_id for this record"""

    fields  = attrb.split(';')
    gene_id = ''
    id_     = ''

    if (len(fields) == 1):
        if ("gene_id" in fields[0]):
            gene_id = fields[0].replace("gene_id", '')
            gene_id = gene_id.strip(' "\n')
        elif ("ID" in fields[0]):
            gene_id = fields[0].replace("ID", '')
            gene_id = gene_id.strip(' "\n=')
        else:
            gene_id = fields[0].strip(' "\n')
        return gene_id
    
    for field in fields:
        field = field.strip(' "')
        if ((field.startswith("gene_id") == False) and (field.startswith("ID") == False)):
            continue
        subfields = field.strip(' "\n,=').split(' ')
        if (len(subfields) == 1 and '=' in subfields[0]):
            subfields = field.strip(' "\n,=').split('=')
        record_id = subfields[-1]
        record_id = record_id.strip(' "\n')

        # hold onto this id if no gene_id found
        if (field[0] == 'I'):
            id_     = record_id
        else:
            gene_id = record_id
            break
            
    if (gene_id == ''):
        gene_id = id_ # assume an ID= was found

    return gene_id

def get_color(num: int) -> list[float]:

    colors = [["cinnabar",   0.89, 0.26, 0.20],
              ["maroon",     0.50, 0, 0],
              ["crimson",    0.86, 0.08, 0.24],      
              ["indigo",     0.29, 0, 0.51],
              ["lavender",   0.71, 0.49, 0.86],
              ["mango",      0.99, 0.75, 0.01], 
              ["peach",      1, 0.90, 0.71],
              ["periwinkle", 0.80, 0.80, 1],
              ["sapphire",   0.06, 0.32, 0.73],
              ["cyan",       0,1,1]]

    return colors[num][1:]


def parse_gtf() -> dict[str, FEATURE]:
    """return a map of features"""

    global gtf, leftbound, rightbound, chrom

    if (gtf.endswith(".gz")):
        fh = gzip.open(gtf, "rt")
    else:
        fh = open(gtf, 'r')

    exon_count  = {}
    genes_found = list()
    gene_bounds = {}
    found_gene  = False
    exons       = {}

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.strip('\n').split('\t')
        if (fields[0] == chrom):
            if (fields[2] == "gene"):
                left_pos  = int(fields[3])
                right_pos = int(fields[4])
                if (leftbound <= left_pos <= rightbound):
                    if right_pos <= rightbound:
                        pass
                    else:
                        right_pos = rightbound
                    found_gene             = True
                    gene_name              = get_gene_id(fields[8])
                    exon_count[gene_name]  = 0
                    gene_bounds[gene_name] = [left_pos, right_pos]
                elif (leftbound <= right_pos <= rightbound):
                    left_pos               = leftbound
                    found_gene             = True
                    gene_name              = get_gene_id(fields[8])
                    exon_count[gene_name]  = 0
                    gene_bounds[gene_name] = [left_pos, right_pos]
                else:
                    found_gene = False
            elif (found_gene):
                if (fields[2] == "exon"):
                    left_pos  = int(fields[3])
                    right_pos = int(fields[4])
                    gene_name = get_gene_id(fields[8])
                    gene_bd   = gene_bounds[gene_name]
                    b1        = gene_bd[0]
                    b2        = gene_bd[1]
                    if (right_pos < b1 or left_pos > b2):
                        continue
                    elif (left_pos < b1 and right_pos < b2):
                        left_pos = b1
                    elif (left_pos < b2 and right_pos > b2):
                        right_pos = b2
                    # adjust for the plot
                    left_pos  = left_pos - leftbound
                    right_pos = right_pos - leftbound
                    if (left_pos < 0 or right_pos < 0):
                        continue
                    if (gene_name not in genes_found):
                        genes_found.append(gene_name)
                    num_genes   = len(genes_found)
                    num         = num_genes % 5
                    gene_color  = get_color(num)
                    direction   = fields[6]
                    exon_count[gene_name] += 1
                    ex  = FEATURE(left_pos, right_pos, direction, gene_color)
                    cnt = exon_count[gene_name]
                    exons[f"{gene_name}_E{cnt}"] = ex
    # print gene names to console
    genes = ' '.join(genes_found)
    print("Found the following genes: " + genes)

    return exons

def parse_gap_file() -> list[FEATURE]:
    """parse a klumpy generated gap file"""

    global gap_file

    if (gap_file == ''):
        return list()

    gaps      = list()
    gap_color = [48/255, 48/255, 48/255]

    if (gap_file.endswith(".gz")):
        fh = gzip.open(gap_file, "rt")
    else:
        fh = open(gap_file, 'r')

    for line in fh:
        fields = line.strip('\n').split('\t')
        if (fields[0] == chrom):
            in_region = False
            left_pos  = int(fields[1])
            right_pos = int(fields[2])
            if (leftbound <= left_pos <= rightbound):
                in_region = True
                if right_pos <= rightbound:
                    pass
                else:
                    right_pos = rightbound
            elif (leftbound <= right_pos <= rightbound):
                left_pos  = leftbound
                in_region = True
            if (in_region):
                left_pos  = left_pos - leftbound
                right_pos = right_pos - leftbound
                gap       = FEATURE(left_pos, right_pos, '+', gap_color)
                gaps.append(gap)

    return gaps            

########## helper functions #############

def get_y_spacing(svg_height: int, y_offset: float, max_cmh: float) -> tuple[float, float]:
    """adjust image height for # of seqs"""

    max_cmh = math.ceil(max_cmh)
    ceiling = copy.copy(max_cmh)

    if (ceiling % 5 == 0):
        ceiling += 5
    else:
        while (ceiling % 5 != 0):
            ceiling += 1

    num_ticks = (ceiling // 5) + 1 # include tick for 0
    space     = svg_height - y_offset
    interval  = space / num_ticks

    return interval, ceiling

########## plotting functions #############


def plot_chrom(cairo_context, seq_length: int, seq_id: str, allocated_pos: float, seq_height: float, x_pos: float):
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
        "Helvetica", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL
    )
    cairo_context.move_to(x, y)
    left = "{:,}".format(leftbound)
    right = "{:,}".format(rightbound)
    seq_id = f"{sequence_id}: {left}-{right}"
    cairo_context.show_text(seq_id)
    cairo_context.stroke()


def create_ticks(seq_length: int) -> tuple[int, int]:
    """Find the appropriate marker lengths for tick marks"""

    marker_sets = {
        1e3:   ["bp", 100],
        5e3:   ["K", 1000],
        7.5e3: ["K", 1500],
        1e4:   ["K", 2000],
        2.5e4: ["K", 5000],
        5e4:   ["K", 10000],
        9e4:   ["K", 15000],
        1e5:   ["K", 20000],
        5e5:   ["K", 25000],
        9e5:   ["K", 30000],
        1e6:   ["M", 50000],
        5e6:   ["M", 70000],
        9e6:   ["M", 100000],
        1e7:   ["M", 150000],
        3e7:   ["M", 200000],
        5e7:   ["M", 500000],
        7e7:   ["M", 700000],
        9e7:   ["M", 900000],
        1e8:   ["M", 1000000],
        5e8:   ["M", 1500000],
        9e8:   ["M", 3000000],
        1e9:   ["M", 5000000],
        }

    multiplier = None

    for marker_length, fields in marker_sets.items():
        if (seq_length <= marker_length):
            multiplier = fields[1]
            break

    # multiplier will be used to show increments in bp length
    if (multiplier != None):
        tick_ct = math.floor(seq_length/multiplier)
    else:
        sys.exit(f"Software was not prepared for sequence/segment of length {seq_length}")

    return tick_ct, multiplier


def draw_ticks(cairo_context, y_coordinate: float, seq_length: int, image_width: float, 
               section: float, x_pos: float) -> int:
    """add tick marks to the reference sequence"""

    global leftbound

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
        if (marker_num < 1e3):
            marker_num = marker_num/100
            marker = "bp"
        elif (1e3 < marker_num < 1e6):
            marker_num = marker_num/1e3
            marker = "K"
        elif (marker_num > 1e6):
            marker_num = marker_num/1e6
            marker = "M"
        marker_num   = round(marker_num, 2)
        marker_label = f"{str(marker_num)}{marker}"
        marker_pos   = x_pos + (allocated_pos * i)
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

    return 0


def plot_label(cairo_context, y: float, theta: float, label: str, midpoint: float) -> int:
    """plot the text labels"""

    new_x = midpoint
    new_y = y - 90 # 63
    cairo_context.save()
    # create labels
    cairo_context.select_font_face("Helvetica", cairo.FONT_SLANT_NORMAL,
                                   cairo.FONT_WEIGHT_NORMAL)
    cairo_context.set_source_rgb(0, 0, 0)  # black
    cairo_context.set_font_size(25)
    cairo_context.move_to(new_x, new_y)
    cairo_context.rotate(theta)
    cairo_context.show_text(label)
    cairo_context.stroke()
    cairo_context.restore()

    return 0

def draw_y_axis(cairo_context, x_pos: float, svg_height: float, y_offset: float, 
                spacing: float, ceiling: float) -> int:
    """just draws a straight line and horizontal ticks"""

    num_ticks = int(ceiling // 5) + 1
    start_y   = svg_height - (y_offset * 1.75) 
    top_y     = start_y - ((num_ticks - 1) * spacing)

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
            "Helvetica", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL
        )
        num_str = str(n * 5)
        num_dim = cairo_context.text_extents(num_str)
        cairo_context.move_to(x_pos - 20 - num_dim.width, y_pos + (num_dim.height / 2))
        cairo_context.show_text(num_str)
        cairo_context.stroke()

    # now add a y-axis label
    y_mid = (start_y - top_y) * 0.65
    ylab = "-log(Q-value)"
    theta = math.radians(270)
    cairo_context.set_source_rgb(0, 0, 0)
    cairo_context.set_font_size(90)
    cairo_context.select_font_face(
        "Helvetica", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL
    )
    cairo_context.move_to(x_pos - 200, y_mid)
    cairo_context.rotate(theta)
    cairo_context.show_text(ylab)
    cairo_context.stroke()

    return 0

def plot_exons(cairo_context, y_coordinate: float, allocated_pos: float, 
               exons: dict[str, FEATURE], x_pos: float) -> int:
    """plot the exonic features"""

    for ex, exon in exons.items():
        fill_in_matches(cairo_context, y_coordinate, x_pos,
                        allocated_pos, exon.lpos, exon.rpos, exon.color)
        # get direction to see if clump is on forward or reverse
        direction = exon.direction
        if (direction == "+"):
            yLabel = y_coordinate + 60
            angle = -45
        elif (direction == "-"):
            yLabel = y_coordinate + 105
            angle = 45
        # add label
        label = ex
        # convert to radians
        theta = math.radians(angle)
        # get the mid point of the text if i were not rotated
        start_pos = x_pos + (exon.lpos * allocated_pos)
        end_pos   = x_pos + (exon.rpos * allocated_pos)
        midpoint  = (end_pos + start_pos) / 2
        # create labels
        plot_label(cairo_context, yLabel, theta, label, midpoint)

    return 0

def fill_in_matches(cairo_context, y: float, x_pos: float, allocated_pos: float, lpos: float, 
                    rpos: float, color: list[float]) -> int:
    """fill in seq with annotation positions"""

    for i in range(lpos, rpos + 1):
        site = x_pos + (i * allocated_pos)  # scale
        cairo_context.move_to(site, y - 20)
        cairo_context.line_to(site, y + 20)
        cairo_context.set_source_rgb(color[0], color[1], color[2]) 
        cairo_context.set_line_width(5)
        cairo_context.stroke()

    return 0


def plot_cmh(cairo_context, snps: list[SNP], svg_height: float, y_offset: float, 
             x_start: float, allocated_pos: float, ceiling: float) -> int:
    """draw the SNPs, which are circles representing the CMH value"""
    allocated_pos_y = svg_height - (y_offset * 1.75)

    for snp in snps:
        x_pos = x_start + ((snp.pos - leftbound) * allocated_pos)
        y_pos = ((1 - (snp.value / ceiling)) * allocated_pos_y) - 15
        cairo_context.arc(x_pos, y_pos, 12, 0, 2*math.pi)
        cairo_context.close_path()
        cairo_context.set_source_rgb(0, 0, 128)
        cairo_context.fill()

    return 0

def draw_image(exons: dict[str, FEATURE], snps: list[SNP], gaps: list[FEATURE], max_cmh: float) -> int:
    """plot clumps using plotting functions"""

    global leftbound, rightbound, chrom

    # this will create a svg for later tweaking
    svg_width    = 5000  # 7500
    svg_height   = 3000
    y_offset     = 250
    chrom_y_pos  = svg_height - y_offset
    chrom_x_pos  = 450
    y_axix_x_pos = 400
    section      = svg_height * 0.02
    seq_len      = rightbound - leftbound
    
    spacing, ceiling = get_y_spacing(svg_height, y_offset, max_cmh)
    
    # create the drawing surface
    ims           = cairo.SVGSurface(f"{chrom}-{leftbound}-{rightbound}.svg", svg_width, svg_height)
    cairo_context = cairo.Context(ims)
    image_width   = svg_width - (chrom_x_pos * 1.5)  # offset by starting positions on both sides

    # create allocated posititions for sites
    allocated_pos = image_width / seq_len

    # plot seqs
    plot_chrom(cairo_context, seq_len, chrom, allocated_pos, chrom_y_pos, chrom_x_pos)
    draw_ticks(cairo_context, chrom_y_pos, seq_len, image_width, section, chrom_x_pos)
    plot_exons(cairo_context, chrom_y_pos, allocated_pos, exons, chrom_x_pos)
    if (len(gaps) > 0):
        for gap in gaps:
            fill_in_matches(cairo_context, chrom_y_pos, chrom_x_pos, allocated_pos, gap.lpos, gap.rpos, gap.color)
    plot_cmh(cairo_context, snps, svg_height, y_offset, chrom_x_pos, allocated_pos, ceiling)
    draw_y_axis(cairo_context, y_axix_x_pos, svg_height, y_offset, spacing, ceiling)
    
    # finish plot
    ims.finish()
    ims.flush()

    return 0

def main() -> int:
    """Parse, plot, and enjoy the visual"""

    # get arguments
    get_arguments()

    # parse inputs
    exons         = parse_gtf()
    snps, max_cmh = parse_cmh()
    gaps          = parse_gap_file()

    # run the entire script
    draw_image(exons, snps, gaps, max_cmh)

    return 0

if __name__ == "__main__":
    main()
