import math, random, sys

from information.predefined import charges,molmass,solventParticles
from utils.abbreviation import R
from utils.parse_command import parse_mol


def solvent(variables,options,lipids,pbc,rest_args):
    pbcx, pbcy, pbcz = pbc

    ################
    ## 3. SOLVENT ##
    ################

    # Charge of the system so far

    last = None
    mcharge = 0
    for j in rest_args["membrane"].atoms:
        if not j[0].strip().startswith('v') and j[1:3] != last:
            mcharge += charges.get(j[1].strip(), 0)
        last = j[1:3]

    last = None
    pcharge = 0
    for j in rest_args['protein'].atoms:
        if not j[0].strip().startswith('v') and j[1:3] != last:
            pcharge += charges.get(j[1].strip(), 0)
        last = j[1:3]

    mcharge = sum([charges.get(i[0].strip(), 0) for i in set([j[1:3] for j in rest_args["membrane"].atoms])])
    pcharge = sum([charges.get(i[0].strip(), 0) for i in
                   set([j[1:3] for j in rest_args['protein'].atoms if not j[0].strip().startswith('v')])])

    charge = mcharge + pcharge
    plen, mlen, slen = 0, 0, 0
    plen = rest_args['protein'] and len(rest_args['protein']) or 0
    mlen = rest_args["membrane"] and len(rest_args["membrane"]) or 0

    print("; NDX Solute %d %d" % (1, rest_args['protein'] and plen or 0))
    print("; Charge of protein: %f" % pcharge)

    print("; NDX Membrane %d %d" % (1 + plen, rest_args["membrane"] and plen + mlen or 0))
    print("; Charge of membrane: %f" % mcharge)
    print("; Total charge: %f" % charge)

    if variables['solv']:

        # First get names and relative numbers for each solvent
        solnames, solnums = zip(*[parse_mol(i) for i in variables['solv']])
        solnames, solnums = list(solnames), list(solnums)
        totS = float(sum(solnums))

        gasnames, gasnums = zip(*[parse_mol(i) for i in variables['gas']])
        gasnames, gasnums = list(gasnames), list(gasnums)
        totG = float(sum(gasnums))

        # Set up a grid
        d = 1 / options["-sold"].value
        dds = options["-sold"].value

        nx, ny, nz = int(1 + d * pbcx), int(1 + d * pbcy), int(1 + d * pbcz)
        dx, dy, dz = pbcx / nx, pbcy / ny, pbcz / nz
        excl, hz = int(nz * options["-excl"].value / pbcz), int(0.5 * nz)

        zshift = 0

        if rest_args["membrane"]:
            memz = [i[2] for i in rest_args["membrane"].coord]
            midz = (max(memz) + min(memz)) / 2
            hz = int(nz * midz / pbcz)  # Grid layer in which the membrane is located
            zshift = (hz + 0.5) * nz - midz  # Shift of membrane middle to center of grid layer

        # Initialize a grid of solvent, spanning the whole cell
        # Exclude all cells within specified distance from membrane center
        grid = [[[i for i in range(nz)] for j in range(ny)] for i in range(nx)]

        if variables['gas']:
            NA = 6.023 * (10 ** 23)
            # Set up a grid
            if not options["-gasd"].value and not options["-gden"].value:
                print("Warning: You have to put in the gas density or the atomic distance")
                sys.exit()

            if options["-gasd"].value:
                fddg = options["-gasd"].value

            if options["-gden"].value:
                density = options["-gden"].value
                for x in gasnames:
                    fddg = pow((molmass[x] / (density * NA * (10 ** (-24)))), (1 / 3))

            if options["-gasd"].value and options["-gden"].value:
                print("Warning: You have entered two quantities of the same meaning")
                sys.exit()

            nxg, nyg, nzg = int(1 + (1 / fddg) * pbcx), int(1 + (1 / fddg) * pbcy), int(1 + (1 / fddg) * pbcz)
            dxg, dyg, dzg = pbcx / nxg, pbcy / nyg, pbcz / nzg
            hzg = int(0.5 * nzg)

            zshiftg = 0

            if rest_args["membrane"]:
                memz = [i[2] for i in rest_args["membrane"].coord]
                midz = (max(memz) + min(memz)) / 2
                hzg = int(nzg * midz / pbcz)  # Grid layer in which the membrane is located
                zshiftg = (hzg + 0.5) * nzg - midz  # Shift of membrane middle to center of grid layer

            # Initialize a grid of solvent, spanning the whole cell
            # Exclude all cells within specified distance from membrane center
            gridg = [[[i for i in range(nzg)] for j in range(nyg)] for i in range(nxg)]

        # Flag all cells occupied by protein or membrane
        for x, y, z in rest_args['protein'].coord + rest_args["membrane"].coord:
            if z >= pbcz:
                x -= rest_args['box'][2][0]
                y -= rest_args['box'][2][1]
                z -= rest_args['box'][2][2]
            if z < 0:
                x += rest_args['box'][2][0]
                y += rest_args['box'][2][1]
                z += rest_args['box'][2][2]
            if y >= pbcy:
                x -= rest_args['box'][1][0]
                y -= rest_args['box'][1][1]
            if y < 0:
                x += rest_args['box'][1][0]
                y += rest_args['box'][1][1]
            if x >= pbcx:
                x -= rest_args['box'][0][0]
            if x < 0:
                x += rest_args['box'][0][0]
            grid[int(nx * x / rest_args['rx'])][int(ny * y / rest_args['ry'])][int(nz * z / rest_args['rz'])] = False
            gridg[int(nxg * x / rest_args['rx'])][int(nyg * y / rest_args['ry'])][int(nzg * z / rest_args['rz'])] = False

        # Flag all cells inside the membrane and put gas inside
        rad = rest_args['r']
        outrad = max(rest_args['maxz']) + rad
        # grid1=[]

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if ((i * dds - pbcx / 2) ** 2 + (j * dds - pbcy / 2) ** 2 + (k * dds - pbcz / 2) ** 2) <= (
                            outrad + 0.5 * excl) ** 2:
                        grid[i][j][k] = False

        for i in range(nxg):
            for j in range(nyg):
                for k in range(nzg):
                    if ((i * fddg - pbcx / 2) ** 2 + (j * fddg - pbcy / 2) ** 2 + (
                            k * fddg - pbcz / 2) ** 2) >= rad ** 2:
                        gridg[i][j][k] = False

        # Set the center for each solvent molecule
        kick = options["-solr"].value
        grid = [(R(), (i + 0.5 + R() * kick) * dx, (j + 0.5 + R() * kick) * dy, (k + 0.5 + R() * kick) * dz)
                for i in range(nx) for j in range(ny) for k in range(nz) if grid[i][j][k]]
        kickg = options["-gasr"].value
        gridg = [(R(), (i + 0.5 + R() * kickg) * dxg, (j + 0.5 + R() * kickg) * dyg, (k + 0.5 + R() * kickg) * dzg)
                 for i in range(nxg) for j in range(nyg) for k in range(nzg) if gridg[i][j][k]]

        # Sort on the random number
        grid.sort()
        gridg.sort()

        # 'grid' contains all positions on which a solvent molecule can be placed.
        # The number of positions is taken as the basis for determining the salt concentration.
        # This is fine for simple salt solutions, but may not be optimal for complex mixtures
        # (like when mixing a 1M solution of this with a 1M solution of that

        # Set the number of ions to add
        nna, ncl = 0, 0
        if options["-salt"]:

            # If the concentration is set negative, set the charge to zero
            if options["-salt"].value.startswith("-"):
                charge = 0
                options["-salt"].value = -float(options["-salt"].value)
            else:
                options["-salt"].value = float(options["-salt"].value)

            # Determine charge to use, either determined or given on command line
            if options["-charge"].value != "0":
                charge = (options["-charge"].value != "auto") and int(options["-charge"].value) or charge
            else:
                charge = 0

            # Determine number of sodium and chloride to add
            concentration = options["-salt"].value
            nsol = ("SPC" in solnames and 1 or 4) * len(grid)
            ncl = max(max(0, charge), int(.5 + .5 * (concentration * nsol / (27.7 + concentration) + charge)))
            nna = ncl - charge

        # Correct number of grid cells for placement of solvent
        ngrid = len(grid) - nna - ncl
        num_sol = [int(ngrid * i / totS) for i in solnums]
        num_gas = [int(len(gridg) * i / totG) for i in gasnums]

        # Add salt to solnames and num_sol
        if nna:
            solnames.append("NA+")
            num_sol.append(nna)
            variables['solv'].append("NA+")
        if ncl:
            solnames.append("CL-")
            num_sol.append(ncl)
            variables['solv'].append("CL-")

        # Names and grid positions for solvent molecules
        solvent = list(zip([s for i, s in zip(num_sol, solnames) for j in range(i)], grid))
        gasin = list(zip([s for i, s in zip(num_gas, gasnames) for j in range(i)], gridg))

        # Extend the list of molecules (for the topology)
        rest_args['molecules'].extend(list(zip(solnames, num_sol)))
        rest_args['molecules'].extend(list(zip(gasnames, num_gas)))

        # Build the solvent
        sol = []
        gas = []
        for resn, (rndm, x, y, z) in solvent:
            rest_args['resi'] += 1
            solmol = solventParticles.get(resn)
            if solmol and len(solmol) > 1:
                # Random rotation (quaternion)
                u, v, w = random.random(), 2 * math.pi * random.random(), 2 * math.pi * random.random()
                s, t = math.sqrt(1 - u), math.sqrt(u)
                qw, qx, qy, qz = s * math.sin(v), s * math.cos(v), t * math.sin(w), t * math.cos(w)
                qq = qw * qw - qx * qx - qy * qy - qz * qz
                for atnm, (px, py, pz) in solmol:
                    qp = 2 * (qx * px + qy * py + qz * pz)
                    rx = x + qp * qx + qq * px + qw * (qy * pz - qz * py)
                    ry = y + qp * qy + qq * py + qw * (qz * px - qx * pz)
                    rz = z + qp * qz + qq * pz + qw * (qx * py - qy * px)
                    sol.append(("%5d%-5s%5s%5d" % (rest_args['resi'] % 1e5, resn, atnm, rest_args['atid'] % 1e5), (rx, ry, rz)))
                    rest_args['atid'] += 1
            else:
                sol.append(
                    ("%5d%-5s%5s%5d" % (rest_args['resi'] % 1e5, resn, solmol and solmol[0][0] or resn, rest_args['atid'] % 1e5), (x, y, z)))
                rest_args['atid'] += 1

        for resn, (rndm, x, y, z) in gasin:
            rest_args['resi'] += 1
            solmol = solventParticles.get(resn)
            if solmol and len(solmol) > 1:
                # Random rotation (quaternion)
                u, v, w = random.random(), 2 * math.pi * random.random(), 2 * math.pi * random.random()
                s, t = math.sqrt(1 - u), math.sqrt(u)
                qw, qx, qy, qz = s * math.sin(v), s * math.cos(v), t * math.sin(w), t * math.cos(w)
                qq = qw * qw - qx * qx - qy * qy - qz * qz
                for atnm, (px, py, pz) in solmol:
                    qp = 2 * (qx * px + qy * py + qz * pz)
                    rx = x + qp * qx + qq * px + qw * (qy * pz - qz * py)
                    ry = y + qp * qy + qq * py + qw * (qz * px - qx * pz)
                    rz = z + qp * qz + qq * pz + qw * (qx * py - qy * px)
                    sol.append(("%5d%-5s%5s%5d" % (rest_args['resi'] % 1e5, resn, atnm, rest_args['atid'] % 1e5), (rx, ry, rz)))
                    rest_args['atid'] += 1
            else:
                sol.append(
                    ("%5d%-5s%5s%5d" % (rest_args['resi'] % 1e5, resn, solmol and solmol[0][0] or resn, rest_args['atid'] % 1e5), (x, y, z)))
                rest_args['atid'] += 1

    rest_args['solvent'] = solvent
    rest_args['plen'] = plen
    rest_args['mlen'] = mlen
    rest_args['sol'] = sol


    return variables,options,lipids,pbc,rest_args
