""" creating venn diagrams"""


from itertools import chain
from collections.abc import Iterable
from matplotlib.patches import Circle, Ellipse
import pylab
# --------------------------------------------------------------------
alignment = {'horizontalalignment': 'center', 'verticalalignment': 'baseline'}


# -------------------------------------------------------------------

def venn(data, outfile, names=None, fill="number", **kwds):
    """
    data: a list
    names: names of groups in data
    fill = ["number"|"logic"|"both"], fill with number, logic label, or both
    outfile = [True|False]
     = [True|False]
    """

    if data is None:
        raise Exception("No data!")
    if len(data) == 2:
        venn2(data, outfile, names, fill,**kwds)
    elif len(data) == 3:
        venn3(data, outfile, names, fill,**kwds)
    elif len(data) == 4:
        venn4(data, outfile, names, fill,**kwds)
    else:
        raise Exception("currently only 2-4 sets venn diagrams are supported")


# --------------------------------------------------------------------
def get_labels(data, fill="number"):
    """
    to get a dict of labels for groups in data
    input
      data: data to get label for
      fill = ["number"|"logic"|"both"], fill with number, logic label, or both
    return
      labels: a dict of labels for different sets
    example:
    In [12]: get_labels([range(10), range(5,15), range(3,8)], fill="both")
    Out[12]:
    {'001': '001: 0',
     '010': '010: 5',
     '011': '011: 0',
     '100': '100: 3',
     '101': '101: 2',
     '110': '110: 2',
     '111': '111: 3'}
    """

    N = len(data)

    sets_data = [set(data[i]) for i in range(N)]  # sets for separate groups
    s_all = set(chain(*data))                             # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value

    labels = {k: "" for k in set_collections}
    if fill == "number":
        for k, val in set_collections.items():
            labels[k] = str(len(val))  # set_collections[k]) for k in set_collections}
    else:  # invalid value
        raise Exception("invalid value for fill")

    return labels

# --------------------------------------------------------------------


def venn2(data, outfile, names=None, fill="number", **kwds):
    """2 circle venn"""

    if (data is None) or len(data) != 2:
        raise Exception("length of data should be 2!")
    if (names is None) or (len(names) != 2):
        names = ("set 1", "set 2")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else:  # default figure size
        figsize = (8, 8)

    fig = pylab.figure(figsize=figsize)
    ax = fig.gca()
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 8)

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 2:
        colors = kwds['colors']
    else:
        colors = ['red', 'green']

    c1 = Circle((3.0, 4.0), radius=2.0, alpha=0.5, color=colors[0])
    c2 = Circle((5.0, 4.0), radius=2.0, alpha=0.5, color=colors[1])

    ax.add_patch(c1)
    ax.add_patch(c2)

    # draw text
    # 1
    pylab.text(2.0, 4.0, labels['10'], **alignment)
    pylab.text(2.5, 4.0, labels['01'], **alignment)
    # 2
    pylab.text(4.0, 4.0, labels['11'], **alignment)
    # names of different groups
    if outfile:
        pylab.text(3.0, 4.0-1.2*2.0, names[0], fontsize=16, **alignment)
        pylab.text(5.0, 4.0-1.2*2.0, names[1], fontsize=16, **alignment)

    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    pylab.savefig(outfile)
# --------------------------------------------------------------------


def venn3(data, outfile, names=None, fill="number", **kwds):
    """3 circle venn"""
    if (data is None) or len(data) != 3:
        raise Exception("length of data should be 3!")
    if (names is None) or (len(names) != 3):
        names = ("set 1", "set 2", "set 3")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else:  # default figure size
        figsize = (10, 10)

    fig = pylab.figure(figsize=figsize)   # set figure size
    ax = fig.gca()
    ax.set_aspect("equal")                # set aspect ratio to 1
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 8)

    # set colors for different Circles or ellipses

    colors = ['red', 'green', 'blue']

    c1 = Circle((3.0, 3.0), radius=2.0, alpha=0.5, color=colors[0])
    c2 = Circle((5.0, 3.0), radius=2.0, alpha=0.5, color=colors[1])
    c3 = Circle((2.5, 4.7), radius=2.0, alpha=0.5, color=colors[2])
    for c in (c1, c2, c3):
        ax.add_patch(c)

    # draw text
    # 1
    pylab.text(3.0-2/2, 3.0-2/2, labels['100'], **alignment)
    pylab.text(5.0+2/2, 3.0-2/2, labels['010'], **alignment)
    pylab.text((3.0+5.0)/2, 4.73+2/2, labels['001'], **alignment)
    # 2
    pylab.text((3.0+5.0)/2, 3.0-2/2, labels['110'], **alignment)
    pylab.text(3.0, 3.0+2*2/3, labels['101'], **alignment)
    pylab.text(5.0, 3.0+2*2/3, labels['011'], **alignment)
    # 3
    pylab.text((3.0+5.0)/2, 3.0+2/3, labels['111'], **alignment)
    # names of different groups
    if outfile:
        pylab.text(1.0, 1.0, names[0], fontsize=16, **alignment)
        pylab.text(7.0, 1.0, names[1], fontsize=16, **alignment)
        pylab.text(4.0, 4.73+1.2*2, names[2], fontsize=16, **alignment)

    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    pylab.savefig(outfile)


def venn4(data, outfile, names=None, fill="number", **kwds):
    """4 circle venn"""

    names = ("set 1", "set 2", "set 3", "set 4")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else:  # default figure size
        figsize = (10, 10)

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 4:
        colors = kwds['colors']
    else:
        colors = ['r', 'g', 'b', 'c']

    # draw ellipse, the coordinates are hard coded in the rest of the function
    fig = pylab.figure(figsize=figsize)   # set figure size
    ax = fig.gca()
    patches = []
    width, height = 170, 110  # width and height of the ellipses
    patches.append(Ellipse((170, 170), width, height, -45, color=colors[0], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -45, color=colors[1], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -135, color=colors[2], alpha=0.5))
    patches.append(Ellipse((230, 170), width, height, -135, color=colors[3], alpha=0.5))
    for e in patches:
        ax.add_patch(e)
    ax.set_xlim(80, 320)
    ax.set_ylim(80, 320)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")

    # draw text
    # 1
    pylab.text(120, 200, labels['1000'], **alignment)
    pylab.text(280, 200, labels['0100'], **alignment)
    pylab.text(155, 250, labels['0010'], **alignment)
    pylab.text(245, 250, labels['0001'], **alignment)
    # 2
    pylab.text(200, 115, labels['1100'], **alignment)
    pylab.text(140, 225, labels['1010'], **alignment)
    pylab.text(145, 155, labels['1001'], **alignment)
    pylab.text(255, 155, labels['0110'], **alignment)
    pylab.text(260, 225, labels['0101'], **alignment)
    pylab.text(200, 240, labels['0011'], **alignment)
    # 3
    pylab.text(235, 205, labels['0111'], **alignment)
    pylab.text(165, 205, labels['1011'], **alignment)
    pylab.text(225, 135, labels['1101'], **alignment)
    pylab.text(175, 135, labels['1110'], **alignment)
    # 4
    pylab.text(200, 175, labels['1111'], **alignment)
    # names of different groups
    if outfile:
        pylab.text(110, 110, names[0], fontsize=16, **alignment)
        pylab.text(290, 110, names[1], fontsize=16, **alignment)
        pylab.text(130, 275, names[2], fontsize=16, **alignment)
        pylab.text(270, 275, names[3], fontsize=16, **alignment)

    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    pylab.savefig(outfile)
