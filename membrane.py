import random,copy
from sympy import symbols

from utils.algebra import *
from utils.structure import Structure
from utils.parse_command import parse_mol

def membrane(variables,options,lipids,pbc,rest_args):
    #################
    ## 2. MEMBRANE ##
    #################
    lipidsx, lipidsy, lipidsz, lipidsa = lipids
    pbcx, pbcy, pbcz = pbc


    membrane = Structure()

    if variables["lipU"]:
        if not options["-r"].value:
            r = pbcx / (2 * math.pi)
        if 2 * options["-r"].value > options["-x"].value:
            r = pbcx / (2 * math.pi)
            print("Warning: the membrane radius is too large ")
        else:
            r = options["-r"].value

        # Create uniformly distributed rectangular coordinates (已改为球面坐标)
        grid = balance(r, rest_args['a'])
        gridcopy = copy.deepcopy(grid)

        # Number of lipids
        print(";lipids number: %d" % len(grid))

        # If there is a protein, mark the corresponding cells as occupied
        if rest_args['protein']:
            for i in gridcopy:
                for j in rest_args['wholecenter']:
                    g = trans1(i, r)
                    center = j[0]
                    newcenter = trans3(trans4(center[0], center[1]), r)

                    d = distance(r, newcenter, g)
                    p = j[1]

                    if d < max(p[0], p[1]):
                        grid.remove(i)

        # Set the XY coordinates
        # To randomize the lipids we add a random number which is used for sorting
        random.seed()
        upperl = []
        upper = []
        for i in grid:
            j = trans2(r, trans1([i], r))
            upperl.append((random.random(), j))
        # Sort on the random number
        upperl.sort()

        for i in upperl:
            upper.append(i[1][0])

        # Upper leaflet
        variables["lipU"], numU = zip(*[parse_mol(i) for i in variables["lipU"]])
        totU = float(sum(numU))
        num_up = [int(len(upper) * i / totU) for i in numU]
        lip_up = [l for i, l in zip(num_up, variables["lipU"]) for j in range(i)]
        leaf_up = (1, zip(lip_up, upper))

        molecules = list(zip(variables["lipU"], num_up))

        leaflet = leaf_up[0]
        leaf_lip = leaf_up[1]

        kick = options["-rand"].value

        coord = []
        maxz = []

        # Build the membrane
        for i in leaf_lip:
            lipid = i[0]
            pos = i[1]
            # Increase the residue number by one
            rest_args['resi'] += 1
            # Set the random rotation for this lipid
            rangle = random.random() * math.pi
            rcos = math.cos(rangle)
            rsin = math.sin(rangle)
            # Fetch the atom list with x,y,z coordinates
            atoms = zip(lipidsa[lipid][1].split(), lipidsx[lipidsa[lipid][0]], lipidsy[lipidsa[lipid][0]],
                        lipidsz[lipidsa[lipid][0]])
            j = []
            for i in atoms:
                if i[0] != "-":
                    j.append(i)
            # Only keep atoms appropriate for the lipid
            at, ax, ay, az = zip(*j)
            # The z-coordinates are spaced at 0.3 nm,
            # starting with the first bead at 0.15 nm
            az = [leaflet * (0.5 + (i - min(az))) * 0.3 for i in az]
            maxz.append(max(az))
            xx = zip(ax, ay)

            nx = []
            ny = []

            newx = []
            newy = []

            cx = []
            cy = []
            cz = []

            for i, j in xx:
                nx.append(rcos * i - rsin * j + pos[0] + random.random() * kick)
                ny.append(rsin * i - rcos * j + pos[1] + random.random() * kick)

            a = pos[0]
            b = pos[1]
            c = trans3(trans4(a, b), r)
            theta = c[1]
            k = 1 / (math.sqrt(2) * pow(math.sin(theta / 2), 2))

            for i in zip(nx, ny):
                x = i[0]
                y = i[1]

                x1 = symbols("x1")
                y1 = symbols("y1")

                x1 = k * x + (1 - k) * a
                y1 = k * y + (1 - k) * b

                newx.append(x1)
                newy.append(y1)

            for ix, iy, iz in zip(newx, newy, az):
                z = trans4(ix, iy)
                phi, theta = trans3(z, r)

                ix = r * math.sin(theta) * math.cos(phi)
                nz = r * math.cos(theta)
                iy = r * math.sin(theta) * math.sin(phi)

                ix += iz * math.sin(theta) * math.cos(phi)
                nz += iz * math.cos(theta)
                iy += iz * math.sin(theta) * math.sin(phi)

                cx.append(ix)
                cy.append(iy)
                cz.append(nz)

                # Add the atoms to the list
            for i in range(len(at)):
                atom = "%5d%-5s%5s%5d" % (rest_args['resi'], lipid, at[i], rest_args['atid'])
                membrane.coord.append((cx[i], cy[i], cz[i]))
                membrane.atoms.append((at[i], lipid, rest_args['resi'], 0, 0, 0))
                rest_args['atid'] += 1

        sumx = 0
        sumy = 0
        sumz = 0

        for i in range(len(membrane.coord)):
            sumx = sumx + membrane.coord[i][0]
            sumy = sumy + membrane.coord[i][1]
            sumz = sumz + membrane.coord[i][2]

        middle = [sumx / len(membrane.coord), sumy / len(membrane.coord), sumz / len(membrane.coord)]

        d1 = math.sqrt((pbcx / 2 - middle[0]) ** 2)
        d2 = math.sqrt((pbcy / 2 - middle[1]) ** 2)
        d3 = math.sqrt((pbcz / 2 - middle[2]) ** 2)

        membrane += (d1, d2, d3)

        rest_args['membrane'] = membrane
        rest_args['r'] = r
        rest_args['molecules'] = molecules
        rest_args['maxz'] = maxz
        variables['numU'] = numU

        return variables,options,lipids,pbc,rest_args