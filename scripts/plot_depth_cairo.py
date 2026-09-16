#!/bin/python

import math
import cairo
import argparse
import sys
import gzip
import copy
import os

# establish some global variables
gtf        = ''
d1         = ''
d2         = ''
d3         = ''
d4         = ''
chrom      = ''
gap_file   = ''
leftbound  = 0
rightbound = 10_000_00
iformat    = "svg"
max_depth  = 0

def get_arguments() -> int:
    """get the arguments"""

    global gtf, d1, d2, d3, d4, gap_file, chrom, leftbound, rightbound, iformat

    formats = ["svg", "pdf", "png"]

    parser = argparse.ArgumentParser(description="Plot a genomic region with sequence coverage for 4 files of coverage")
    parser.add_argument("-g", "--gtf", help="Annotations file in GTF/GFF3 format", required=True)
    parser.add_argument("-1", "--d1", help="depth file", required=True)
    parser.add_argument("-2", "--d2", help="depth file", required=True)
    parser.add_argument("-3", "--d3", help="depth file", required=True)
    parser.add_argument("-4", "--d4", help="depth file", required=True)
    parser.add_argument("-G", "--gap_file", help="gaps file generated via Klumpy[optional]", default='')
    parser.add_argument("-l", "--leftbound", help="lefbound", default=leftbound, type=int)
    parser.add_argument("-r", "--rightbound", help="rightbound", default=rightbound,type=int)
    parser.add_argument("-C", "--chrom", help="name of chromosome to plot", default=chrom, required=True)
    parser.add_argument("-f", "--format", help="output image format [default: svg]", default=iformat, choices=formats)

    args       = parser.parse_args()
    gtf        = args.gtf
    d1         = args.d1
    d2         = args.d2
    d3         = args.d3
    d4         = args.d4
    gap_file   = args.gap_file
    leftbound  = args.leftbound
    rightbound = args.rightbound
    chrom      = args.chrom
    iformat    = args.format

    assert leftbound > 0, "Leftbound must be greater than 0"
    assert leftbound < rightbound, "Leftbound must be smaller than rightbound"
    assert chrom != '', "Chromosome name needs to be specified"

    for f in [gtf, d1, d2, d3, d4]:
        assert os.path.isfile(f), f"Could not locate {f}"

    if (gap_file != ''):
        assert os.path.isfile(gap_file), f"Could not locate {gap_file}"

    return 0

class SITE:
    def __init__(self, pos: int, value: int):
        self.pos   = pos
        self.value = value

class FEATURE:
    def __init__(self, left_pos: int, right_pos: int, direction: int, color: list[float]):
        self.lpos      = left_pos
        self.rpos      = right_pos
        self.direction = direction
        self.color     = color

#################
### functions ###
#################

######### parsing functions ############
def parse_depth(dfile: str) -> list[SITE]:
    """parse a simple tsv file containing sequencing coverage info"""

    global leftbound, rightbound, chrom, max_depth

    positions = list()
    prev      = leftbound
    local_max = 0

    if (dfile.endswith(".gz")):
        fh = gzip.open(dfile, "rt")
    else:
        fh = open(dfile, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue

        fields = line.strip('\n').split('\t')

        if (fields[0] != chrom):
            continue
        pos = int(fields[1])
        if (leftbound <= pos <= rightbound):
            value     = int(fields[2])
            max_depth = max(max_depth, value)
            local_max = max(local_max, value)
            if (prev != pos - 1 and pos != leftbound):
                for i in range(prev, pos):
                    positions.append(SITE(i, 0))
            prev = pos
            positions.append(SITE(pos, value))
        # we assume that the depth file is sorted
        if (pos > rightbound):
            break
    fh.close()

    # print(f"Local max for {os.path.basename(dfile)}: {local_max}")

    return positions

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

    # place it here in main since we actually processed this file
    print(f"Finished parsing {gap_file}")

    return gaps            

########## helper functions #############

def get_y_spacing(img_height: int, y_offset: float) -> tuple[float, float, float]:
    """adjust pdf height for # of seqs"""

    global max_depth

    if (max_depth <= 100):
        divisor = 10 # <-- this was hard coded in
    elif (max_depth <= 200):
        divisor = 50
    else:
        divisor = 100

    max_d    = math.ceil(max_depth)
    ceiling  = copy.copy(max_depth)

    if (ceiling % divisor == 0):
        ceiling += divisor # may want to change this later
    # else:
    #     while (ceiling % divisor != 0):
    #         ceiling += divisor

    num_ticks       = (ceiling // divisor) # include tick for 0
    space           = img_height - y_offset
    interval        = space / num_ticks
    allocated_space = space / max_d

    return interval, ceiling, allocated_space

########## plotting functions #############


def plot_chrom(cairo_context, seq_length: int, seq_id: str, allocated_pos: float, 
               seq_height: float, x_pos: float) -> int:
    """plot a rectange to represent the genomic locus"""

    global leftbound, rightbound

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
    left   = "{:,}".format(leftbound)
    right  = "{:,}".format(rightbound)
    seq_id = f"{sequence_id}: {left}-{right}"
    cairo_context.show_text(seq_id)
    cairo_context.stroke()

    return 0

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
        # print(marker_label)
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

def draw_y_axis(cairo_context, x_pos: float, fig_height: float, y_offset: float, 
                spacing: float, ceiling: float) -> int:
    """just draws a straight line and horizontal ticks"""

    num_ticks = (ceiling // 10)
    start_y   = fig_height - (y_offset * 1.75)
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
        if ((n * 100) < 100):
            xoffset = 100
        elif ((n * 100) >= 1000):
            xoffset = 200
        else:
            xoffset = 150
        cairo_context.move_to(x_pos - xoffset, y_pos)
        cairo_context.show_text(str(n * 100))
        cairo_context.stroke()

    # now add a y-axis label
    ylab = "Read Coverage"
    theta = math.radians(270)
    cairo_context.set_source_rgb(0, 0, 0)
    cairo_context.set_font_size(90)
    cairo_context.select_font_face(
        "Helvetica", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL
    )
    label_dim = cairo_context.text_extents(ylab)
    # y_mid   = (start_y - top_y) * 0.65 
    y_mid     = ((start_y + top_y) / 2) + (label_dim.width / 2)
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

def plot_depth(cairo_context, sites: list[SITE], fig_height: float, y_offset: float, 
               x_start: float, allocated_pos: float, allocated_pos_y: float, 
               ceiling: float, color: list[float]) -> int:
    """draw a single line with each site being a base pair position"""

    global leftbound

    allocated_pos_y = fig_height - (y_offset * 1.75)
    # allocated_pos_y = allocated_pos_y - (((num_ticks - 1) * spacing))

    for i, site in enumerate(sites):
        x_pos = x_start + ((site.pos - leftbound) * allocated_pos)
        y_pos = (1 - (site.value / ceiling)) * allocated_pos_y
        if (i == 0):
            cairo_context.move_to(x_pos, y_pos)
        else:
            cairo_context.line_to(x_pos, y_pos)
    
    cairo_context.set_source_rgba(color[0], color[1], color[2], color[3]) 
    cairo_context.set_line_width(5)
    cairo_context.stroke()

    return 0

def draw_image(exons: dict[str, FEATURE], all_sites: list[list[SITE]], gaps: list[FEATURE]) -> int:
    """the function is the orchestrator in this script"""

    global chrom, leftbound, rightbound, iformat

    fig_width    = 5000  # 7500
    fig_height   = 3000
    y_offset     = 250
    chrom_y_pos  = fig_height - y_offset
    chrom_x_pos  = 450
    y_axix_x_pos = 400
    section      = fig_height * 0.02
    seq_len      = rightbound - leftbound
        
    spacing, ceiling, allocated_space_y = get_y_spacing(fig_height, y_offset)

     # create the drawing surface
    outname = f"{chrom}-{leftbound}-{rightbound}_depth.{iformat}"
    if (iformat == "svg"):
        ims = cairo.SVGSurface(outname, int(fig_width), int(fig_height))
    elif (iformat == "pdf"):
        ims = cairo.PDFMetadata(outname, int(fig_width), int(fig_height))
    elif (iformat == "png"):
        ims = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(fig_width), int(fig_height))

    cairo_context = cairo.Context(ims)
    image_width   = fig_width - (chrom_x_pos * 1.5)  # offset by starting positions on both sides
    allocated_pos = image_width / seq_len # create allocated posititions for sites
    
    # plot the genomic locus & the features
    plot_chrom(cairo_context, seq_len, chrom, allocated_pos, chrom_y_pos, chrom_x_pos)
    draw_ticks(cairo_context, chrom_y_pos, seq_len, image_width, section, chrom_x_pos)
    plot_exons(cairo_context, chrom_y_pos, allocated_pos, exons, chrom_x_pos)
    if (len(gaps) > 0):
        for gap in gaps:
            fill_in_matches(cairo_context, chrom_y_pos, chrom_x_pos, allocated_pos, gap.lpos, gap.rpos, gap.color)
    for i, sites in enumerate(all_sites):
        if (i < 2):
            color = [227/255, 30/255, 51/255, 0.75] # reds will be red-ish
        else:
            color = [79/255, 180/255, 230/255, 0.75] # yellows will be blue-ish
        plot_depth(cairo_context, sites, fig_height, y_offset, chrom_x_pos, allocated_pos, allocated_space_y, ceiling, color)

    draw_y_axis(cairo_context, y_axix_x_pos, fig_height, y_offset, spacing, ceiling)
    
    # finish plot
    if (iformat == "svg"):
        ims.finish()
        ims.flush()
    elif (iformat == "pdf"):
        cairo_context.show_page()
    elif (iformat == "png"):
        ims.write_to_png(outname)

    return 0

def main():
    """Entry point to this application"""

    # get arguments
    get_arguments()

    # parse inputs
    exons = parse_gtf()
    print(f"Finished parsing {gtf}")

    all_sites = []
    for d in [d1, d2, d3, d4]:
        sites = parse_depth(d)
        print(f"Finished parsing {d}")
        all_sites.append(sites)

    gaps = parse_gap_file()
    # run the entire script
    draw_image(exons, all_sites, gaps)

if __name__ == "__main__":
    main()

