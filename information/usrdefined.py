import sys
from .predefined import lipidsx,lipidsy,lipidsz,lipidsa

def user_define(variables,options):

    # HII edit - lipid definition
    # Add specified lipid definition to insane lipid library
    for name, head, link, tail in zip(variables['usrmols'], variables['usrheads'], variables['usrlinks'], variables['usrtails']):
        print("Lipid name %s : %s - %s - %s" % (name, head, link, tail))

        moltype = "usr_" + name
        lipidsx[moltype] = []
        lipidsy[moltype] = []
        lipidsz[moltype] = []
        headArray = (head).split()
        linkArray = (link).split()
        tailsArray = (tail).split()
        lipidDefString = ""

        if len(tailsArray) != len(linkArray):
            print("Error, Number of tails has to equal number of linkers")
            sys.exit()

        # Find longest tail
        maxTail = 0
        for cTail in tailsArray:
            if len(cTail) > maxTail:
                maxTail = len(cTail)
        cBeadZ = maxTail + len(
            headArray)  # longest tail + linker (always x1) + lengths of all heads - 1 (as it starts on 0)

        # Add head beads
        for cHead in headArray:
            lipidsx[moltype].append(0)
            lipidsy[moltype].append(0)
            lipidsz[moltype].append(cBeadZ)
            cBeadZ -= 1
            lipidDefString += variables['usrLipHeadMapp'][cHead] + " "

        # Add linkers
        for i, cLinker in enumerate(linkArray):
            lipidsx[moltype].append(max(i - 0.5, 0))
            lipidsy[moltype].append(0)
            lipidsz[moltype].append(cBeadZ)
            lipidDefString += "GL" + str(i + 1) + " "

        # Add tails
        for i, cTail in enumerate(tailsArray):
            cBeadZ = maxTail - 1

            for j, cTailBead in enumerate(cTail):
                lipidsx[moltype].append(i)
                lipidsy[moltype].append(0)
                lipidsz[moltype].append(cBeadZ)
                cBeadZ -= 1
                lipidDefString += cTailBead + str(j + 1) + variables['usrIndexToLetter'][i] + " "

        lipidsa[name] = (moltype, lipidDefString)
    # End user lipid definition


    # HII edit - lipid definition, had to move this one below the user lipid definitions to scale them to.
    # First all X/Y coordinates of templates are centered and scaled (magic numbers!)
    for i in lipidsx.keys():
        cx = (min(lipidsx[i]) + max(lipidsx[i])) / 2
        lipidsx[i] = [0.25 * (j - cx) for j in lipidsx[i]]
        cy = (min(lipidsy[i]) + max(lipidsy[i])) / 2
        lipidsy[i] = [0.25 * (j - cy) for j in lipidsy[i]]

    # Periodic boundary conditions

    # option -pbc keep really overrides everything
    if options["-pbc"].value == "keep" and variables["tm"]:
        options["-x"].value = variables["tm"][0].box[:3]
        options["-y"].value = variables["tm"][0].box[3:6]
        options["-z"].value = variables["tm"][0].box[6:]

    # options -x, -y, -z take precedence over automatic determination
    pbcSetX = 0
    if type(options["-x"].value) in (list, tuple):
        pbcSetX = options["-x"].value
    elif options["-x"].value:
        pbcSetX = [options["-x"].value, 0, 0]

    pbcSetY = 0
    if type(options["-y"].value) in (list, tuple):
        pbcSetY = options["-y"].value
    elif options["-y"].value:
        pbcSetY = [0, options["-y"].value, 0]

    pbcSetZ = 0
    if type(options["-z"].value) in (list, tuple):
        pbcSetZ = options["-z"].value
    elif options["-z"].value:
        pbcSetZ = [0, 0, options["-z"].value]

    a = options["-a"].value

    return variables,options,(lipidsx,lipidsy,lipidsz,lipidsa),(pbcSetX,pbcSetY,pbcSetZ),{'a':a}