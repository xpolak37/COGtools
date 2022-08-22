import pkg_resources
from PIL import Image, ImageDraw, ImageFont, ImageColor
import os
import pandas as pd


def get_track_template(pos_track=(0.95, 0.90, 0.85, 0.80), size=10.0, cog_palette=True, output_dir=os.getcwd()):
    """
    Generate file for Track Manager option in DNAPlotter
    :type cog_palette: bool
    :param pos_track: the positions for plotting features (CDS forward strand, CDS reverse strand, pseudogenes, RNA genes)
    :param size: the size of track
    :param cog_palette: use palette from COG database
    :param output_dir: the output file
    :return: track template file
    """
    # table for new data
    track_template = pd.DataFrame(columns=["#pos", "size", "for", "rev", "not", "any", "key", "qual", "val", "col"])
    file = pkg_resources.resource_filename(__name__, 'COGtools-data/fun-20.tab.txt')
    features = (open(file, "r")).read()
    features = features.split('\n')
    index = 1 if cog_palette else 3
    for feature in features:
        feature = feature.split('\t')
        if "RNA" not in feature[0] and "-" not in feature[0]:
            color = ImageColor.getcolor('#' + feature[index], "RGB")
            color = str(color[0]) + ":" + str(color[1]) + ":" + str(color[2])
            new_row = pd.DataFrame(
                {"#pos": [pos_track[0],pos_track[1],pos_track[2]],
                 "size": [str(size)] * 3,
                 "for": ["true","false","true"],
                 "rev": ["false","true","true"],
                 "not": ["false"]*3,
                 "any": ["false"]*3,
                 "key": ["CDS","CDS", "pseudogene"],
                 "qual": ["CAT"]*3,
                 "val": [feature[0]]*3,
                 "col": [color]*3})
            track_template = pd.concat([track_template, new_row], ignore_index=True)
        else:
            color = ImageColor.getcolor('#' + feature[1], "RGB")
            color = str(color[0]) + ":" + str(color[1]) + ":" + str(color[2])
            new_row = pd.DataFrame(
                {"#pos": [pos_track[3]],
                 "size": [str(size)],
                 "for": ["true"],
                 "rev": ["true"],
                 "not": ["false"],
                 "any": ["false"],
                 "key": [feature[0]],
                 "qual": ["null"],
                 "val": ["null"],
                 "col": [color]})
            track_template = pd.concat([track_template, new_row], ignore_index=True)

    f = open(output_dir + '/track_template', 'w')
    f.write('# created with COGtools 1.0.0\n'
            '# Track template for DNAPlotter\n')
    f.close()
    return track_template.to_csv(output_dir + '/track_template', sep='\t', index=False, mode='a')


def get_legend(font='arial.ttf',output_dir=os.getcwd(), cog_palette=True):
    """
    Create a legend for the genome map
    """
    # import description for the legend
    file = pkg_resources.resource_filename(__name__, 'COGtools-data/fun-20.tab.txt')
    features = (open(file, "r")).read()
    features = features.split('\n')

    # create a white image
    img = Image.new('RGB', (1700, 2500), 'white')
    image_edit = ImageDraw.Draw(img)
    myFont = ImageFont.truetype(font, 50)
    start = 50

    # add individual objects to the legend
    index = 1 if cog_palette else 3
    for feature in features:
        feature = feature.split('\t')
        text = feature[0] if len(feature) == 2 else feature[0] + '  ' + feature[2]
        if len(feature) == 2:
            color = ImageColor.getcolor('#' + feature[1], "RGB")
        else:
            color = ImageColor.getcolor('#' + feature[index], "RGB")
        text = 'COG unknown' if text == '-' else text
        image_edit.rectangle((50, start, 50+80, start+80), fill=color)
        image_edit.text((150, start+15), text, font=myFont, fill=(0, 0, 0))
        start = start + 80

    img.save(output_dir + '/legend.jpg')
