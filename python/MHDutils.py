import os
from scipy.optimize import fsolve

def tagTag(dict):
    """
    Adds '/' before and after the strings working as keys in a given dictionary
    and returns the modified dictionary
    """
    tags = list(dict.keys())
    values = list(dict.values())
    new_dict = {}
    for i in range(len(tags)):
        tags[i] = '/' + tags[i] + '/'
        new_dict[tags[i]] = values[i]
    else:
        return new_dict


def replaceDictList(tag_dict, filePath, verbose=False):
    """
    Replaces tags in a dictionary 'tag_dict' (tags=keys) present in a given
    file 'filePath' with a list of values in the same dictionary
    """
    tag_dict = tagTag(tag_dict)
    with open(filePath, 'r') as (f):
        fileText = f.read()
    for tag in tag_dict:
        if tag in fileText:
            if verbose:
                print('Replacing {} tag in {} file. Value: {}.'
                      .format(tag, filePath, tag_dict[tag]))
            fileText = fileText.replace(tag, str(tag_dict[tag]))
    else:
        with open(filePath, 'w') as (f):
            f.write(fileText)


def replaceTags(tag_dict, pathMain):
    """
    Replaces all tags using the 'replaceDictList' function
    and looping over the whole directory
    Designed to give 'pathMain' the name of an OpenFOAM case directory
    """
    for f1 in os.listdir(pathMain):
        path2 = pathMain + '/' + f1
        if os.path.isdir(path2):
            for f2 in os.listdir(path2):
                path3 = path2 + '/' + f2
                if os.path.isdir(path3):
                    for f3 in os.listdir(path3):
                        path4 = path3 + '/' + f3
                        replaceDictList(tag_dict, path4)

                else:
                    replaceDictList(tag_dict, path3)


def lenC2CN(len, N=25, c2c=1.05):
    """
    Returns grading, and minimun and maximum cell lenghts for a given length,
    number of cells and cell-to-cell ratio following a geometric progression
    """
    G = c2c ** (N - 1)
    cmin = len * (c2c - 1) / (c2c ** N - 1)
    cmax = cmin * G
    return (G, cmin, cmax)


def lenCminC2C(len, cmin, c2c=1.1, N_ini=2, G_ini=2, verbose=False):

    def equations(p):
        N, G = p
        return (c2c ** N - 1 + len * (1 - c2c) / cmin, c2c ** (N - 1) - G)

    N, G = fsolve(equations, (N_ini, G_ini))
    while abs(c2c ** (N - 1) - G) > 0.001:
        N_ini = N_ini * 2
        G_ini = G_ini * 1.5
        N, G = fsolve(equations, (N_ini, G_ini))

    cmax = cmin * G
    if verbose:
        print('N= ' + str(int(N)) + ' (c2c**(N-1)-G)>0.01 = ' + str(c2c ** (N - 1) - G > 0.01))
    return int(N), G, cmax


def getLatestTime(myDir):
    latestTime = '0'
    for filename in os.listdir(myDir):
        if float(filename) > float(latestTime):
            latestTime = filename
    return latestTime
