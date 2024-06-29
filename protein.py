from utils.structure import Structure
from utils.algebra import *
from utils.abbreviation import R
from information.predefined import apolar

def protein(variables,options,lipids,pbcSet,rest_args):
    pbcSetX, pbcSetY, pbcSetZ = pbcSet

    ################
    ## I. PROTEIN ##
    ################

    protein = Structure()
    wholecenter = []
    prot = []
    shift = [0]  # Shift in x direction per protein

    ## A. NO PROTEIN ---
    if not variables["tm"]:  # tm: [Structure, Structure, ...]

        resi = 0
        # Set the box -- If there is a hole, add its radius to the distance
        pbcx = pbcy = pbcz = options["-d"].value + (options["-hole"] and options["-hole"].value or 0)
        if "hexagonal".startswith(options["-pbc"].value):
            # Hexagonal prism -- y derived from x directly
            pbcy = math.sqrt(3) * pbcx / 2
            pbcz = options["-dz"].value or options["-z"].value or options["-d"].value
        elif "optimal".startswith(options["-pbc"].value):
            # Rhombic dodecahedron with hexagonal XY plane
            pbcy = math.sqrt(3) * pbcx / 2
            pbcz = math.sqrt(6) * options["-d"].value / 3
        if "rectangular".startswith(options["-pbc"].value):
            pbcz = options["-dz"].value or options["-z"].value or options["-d"].value

        # Possibly override
        pbcx = pbcSetX and pbcSetX[0] or pbcx
        pbcy = pbcSetY and pbcSetY[1] or pbcy
        pbcz = pbcSetZ and pbcSetZ[2] or pbcz


    ## B. PROTEIN ---
    else:

        for prot in variables["tm"]:

            ## a. NO MEMBRANE --
            if not variables['lipL']:

                # A protein, but don't add lipids... Just solvate the protein
                # Maybe align along principal axes and then build a cell according to PBC

                # Set PBC starting from diameter and adding distance
                if "cubic".startswith(options["-pbc"].value):
                    pbcx = pbcy = pbcz = prot.diam() + options["-d"].value
                elif "rectangular".startswith(options["-pbc"].value):
                    pbcx, pbcy, pbcz = vvadd(vvsub(prot.fun(max), prot.fun(min)), options["-d"].value)
                else:
                    # Rhombic dodecahedron
                    pbcx = pbcy = prot.diam() + options["-d"].value
                    pbcz = math.sqrt(2) * pbcx / 2

                # Possibly override
                pbcx = pbcSetX and pbcSetX[0] or pbcx
                pbcy = pbcSetY and pbcSetY[1] or pbcy
                pbcz = pbcSetZ and pbcSetZ[2] or pbcz

                # Center coordinates in rectangular brick -- Add solvent next
                if len(variables["tm"]) == 1:
                    prot.center((0.5 * pbcx, 0.5 * pbcy, 0.5 * pbcz))

                # Do not set an exclusion range for solvent
                options["-excl"].value = -1


            ## b. PROTEIN AND MEMBRANE --
            else:

                # Have to build a membrane around the protein.
                # So first put the protein in properly.

                # Center the protein and store the shift
                shift = prot.center((0, 0, 0))

                ## 1. Orient with respect to membrane
                # Orient the protein according to the TM region, if requested   # TM: 跨膜区域
                # This doesn't actually work very well...(因为实际上是用的蛋白质的质心及非极性原子来计算的，可能不准确)
                if options["-orient"]:

                    # Grid spacing (nm)
                    d = options["-od"].value  # 网格最小间距
                    pw = options["-op"].value  # 亲疏水性的权重

                    # Determine grid size
                    mx, my, mz = prot.fun(min)
                    rx, ry, rz = prot.fun(lambda x: max(x) - min(x) + 1e-8)

                    # Number of grid cells
                    nx, ny, nz = int(rx / d + 0.5), int(ry / d + 0.5), int(rz / d + 0.5)

                    # Initialize grids
                    atom = [[[0 for i in range(nz + 2)] for j in range(ny + 2)] for k in range(nx + 2)]
                    phobic = [[[0 for i in range(nz + 2)] for j in range(ny + 2)] for k in range(nx + 2)]
                    surface = []
                    for i, (ix, iy, iz) in zip(prot.atoms, prot.coord):
                        if i[1] != "DUM":  # DUM:虚拟或占位的原子，而不是真实的蛋白质结构中的一部分
                            jx, jy, jz = int(nx * (ix - mx) / rx), int(ny * (iy - my) / ry), int(nz * (iz - mz) / rz)
                            atom[jx][jy][jz] += 1
                            phobic[jx][jy][jz] += (i[1].strip() in apolar)  # 同样亲疏水性也是binary表示

                    # Determine average density
                    occupd = sum([bool(k) for i in atom for j in i for k in j])
                    avdens = float(sum([sum(j) for i in atom for j in i])) / occupd

                    # cgofile  = open('density.cgo',"w")
                    # cgofile.write('[\n')
                    for i in range(nx):
                        for j in range(ny):
                            for k in range(nz):
                                if atom[i][j][k] > 0.1 * avdens:
                                    # Check the neighbouring cells; If one of them is not occupied, count cell as surface
                                    if not (atom[i - 1][j][k] and atom[i + 1][j][k] and
                                            atom[i][j - 1][k] and atom[i][j + 1][k] and
                                            atom[i][j][k - 1] and atom[i][j][k + 1]):
                                        sx, sy, sz = mx + rx * (i + 0.5) / nx, my + ry * (j + 0.5) / ny, mz + rz * (
                                                    k + 0.5) / nz
                                        sw = (2.0 * phobic[i][j][k] / atom[i][j][k]) ** pw
                                        surface.append((sx, sy, sz, sw))
                                        # cgofile.write("    7.0, %f, %f, %f, %f,\n"%(10*sx,10*sy,10*sz,0.25*sw))
                    # cgofile.write(']\n')
                    # cgofile.close()

                    sx, sy, sz, w = zip(*surface)
                    W = 1.0 / sum(w)

                    # Weighted center of apolar region; has to go to (0,0,0)
                    sxm, sym, szm = [sum(p) * W for p in
                                     zip(*[(m * i, m * j, m * k) for m, i, j, k in zip(w, sx, sy, sz)])]

                    # Place apolar center at origin
                    prot.center((-sxm, -sym, -szm))
                    sx, sy, sz = zip(*[(i - sxm, j - sym, k - szm) for i, j, k in zip(sx, sy, sz)])

                    # Determine weighted deviations from centers
                    dx, dy, dz = zip(*[(m * i, m * j, m * k) for m, i, j, k in zip(w, sx, sy, sz)])

                    # Covariance matrix for surface
                    xx, yy, zz, xy, yz, zx = [sum(p) * W for p in
                                              zip(*[(i * i, j * j, k * k, i * j, j * k, k * i) for i, j, k in
                                                    zip(dx, dy, dz)])]

                    # PCA: u,v,w are a rotation matrix
                    (ux, uy, uz), (vx, vy, vz), (wx, wy, wz), r = mijn_eigen_sym_3x3(xx, yy, zz, xy, zx, yz)

                    # Rotate the coordinates
                    prot.coord = [(ux * i + uy * j + uz * k, vx * i + vy * j + vz * k, wx * i + wy * j + wz * k) for
                                  i, j, k in prot.coord]

                ## 4. Orient the protein in the xy-plane
                ## i. According to principal axes and unit cell
                if options["-rotate"].value == "princ":

                    x, y, z = zip(*prot.coord)

                    # The rotation matrix in the plane equals the transpose
                    # of the matrix of eigenvectors from the 2x2 covariance
                    # matrix of the positions.
                    # For numerical stability we do
                    # d_i     = x_i - x_0
                    # mean(x) = x_0 + sum(d_i)/N =
                    # var(x)  = sum((d_i - mean(d))**2)/(N-1)
                    xy = ssd(x, y)
                    if xy != 0:
                        xx = ssd(x, x)
                        yy = ssd(y, y)

                        # The eigenvalues are the roots of the 2nd order
                        # characteristic polynomial, with the coefficients
                        # equal to the trace and the determinant of the
                        # matrix.
                        t, d = xx + yy, xx * yy - xy * xy
                        # The two eigenvectors form a 2D rotation matrix
                        # R = ((cos,sin),(-sin,cos)), which means that
                        # the second eigenvector follows directly from
                        # the first. We thus only need to determine one.
                        l1 = t / 2 + math.sqrt(0.25 * t * t - d)

                        ux, uy = l1 - yy, xy
                        lu = math.sqrt(ux * ux + uy * uy)

                        ux /= lu
                        uy /= lu

                        # Finally we rotate the system in the plane by
                        # matrix multiplication with the transpose of
                        # the matrix of eigenvectors
                        prot.coord = [(ux * i + uy * j, ux * j - uy * i, k) for i, j, k in zip(x, y, z)]

                ## ii. Randomly
                elif options["-rotate"].value == "random":
                    ux = math.cos(R() * 2 * math.pi)
                    uy = math.sqrt(1 - ux * ux)
                    prot.coord = [(ux * i + uy * j, ux * j - uy * i, k) for i, j, k in prot.coord]

                ## iii. Specifically
                elif options["-rotate"]:
                    ux = math.cos(float(options["-rotate"].value) * math.pi / 180.)
                    uy = math.sin(float(options["-rotate"].value) * math.pi / 180.)
                    prot.coord = [(ux * i + uy * j, ux * j - uy * i, k) for i, j, k in prot.coord]

                ## 5. Determine the minimum and maximum x and y of the protein
                pmin, pmax = prot.fun(min), prot.fun(max)
                prng = (pmax[0] - pmin[0], pmax[1] - pmin[1], pmax[2] - pmin[2])
                center = (0.5 * (pmin[0] + pmax[0]), 0.5 * (pmin[1] + pmax[1]))

                wholecenter.append((center, prng))

                # Set the z-dimension
                pbcz = pbcSetZ and pbcSetZ[2]
                # If it is not set, set pbcz to the dimension of the protein
                pbcz = pbcz or prng[2]
                pbcz += options["-dz"].value or options["-d"].value or 0

                # At this point we should shift the subsequent proteins such that they end up
                # at the specified distance, in case we have a number of them to do
                # y-shift is always -ycenter
                # x-shift is -xmin+distance+xmax(current)
                xshft, yshft = shift[-1] - pmin[0] + (options["-d"].value or 0), -center[1]
                shift.append(shift[-1] + pmax[0] + (options["-d"].value or 0))

                ## 6. Set box (brick) dimensions
                pbcx = (options["-d"].value or 0) + prng[0]
                if "square".startswith(options["-pbc"].value):
                    pbcy = pbcx
                elif "rectangular".startswith(options["-pbc"].value):
                    pbcy = options["-d"].value + prng[1]
                else:
                    # This goes for a hexagonal cell as well as for the optimal arrangement
                    # The latter is hexagonal in the membrane plane anyway...
                    pbcy = math.cos(math.pi / 6) * pbcx

                ## 7. Adjust PBC for hole
                # If we need to add a hole, we have to scale the system
                # The scaling depends on the type of PBC
                if options["-hole"]:
                    if ("square".startswith(options["-pbc"].value) or
                            "rectangular".startswith(options["-pbc"].value)):
                        scale = 1 + options["-hole"].value / min(pbcx, pbcy)
                    else:
                        area = options["-hole"].value ** 2 / math.cos(math.pi / 6)
                        scale = 1 + area / (pbcx * pbcy)
                    pbcx, pbcy = scale * pbcx, scale * pbcy

                pbcx = pbcSetX and pbcSetX[0] or pbcx
                pbcy = pbcSetY and pbcSetY[1] or pbcy

                ## 2. Shift of protein relative to the membrane center
                if options["-dm"]:
                    if options["-dm"].value < 0:
                        zshift = options["-dm"].value  # - max(zip(*prot.coord)[2])
                    else:
                        zshift = options["-dm"].value  # - min(zip(*prot.coord)[2])
                elif not options["-center"]:
                    zshift = -shift[2]
                else:
                    zshift = 0

                # Now we center the system in the rectangular
                # brick corresponding to the unit cell
                # If -center is given, also center z in plane
                prot += (0.5 * pbcx, 0.5 * pbcy, zshift)

                rescoord = []

            # And we collect the atoms
            protein.atoms.extend(prot.atoms)
            protein.coord.extend(rescoord)

        # Extract the parts of the protein that are in either leaflet
        # 从蛋白质结构中提取出位于任一脂双层（leaflet）中的部分。
        # 在双层膜结构中，通常存在两个叫做"leaflet"的层，一个在上半部分，一个在下半部分，它们之间由蛋白质通道或其他分子分隔开
        prot_up, prot_lo = [], []
        for ix, iy, iz in protein.coord:
            if iz > 0 and iz < 2.4:
                prot_up.append((ix, iy))
            elif iz < 0 and iz > -2.4:
                prot_lo.append((ix, iy))

        # Current residue ID is set to that of the last atom
        resi = protein.atoms[-1][2]

    atid = len(protein) + 1
    molecules = []

    # The box dimensions are now (likely) set.
    # If a protein was given, it is positioned in the center of the
    # rectangular brick.

    # Set the lattice vectors
    if ("rectangular".startswith(options["-pbc"].value) or
            "square".startswith(options["-pbc"].value) or
            "cubic".startswith(options["-pbc"].value)):
        box = [[pbcx, 0, 0], [0, pbcy, 0], [0, 0, pbcz]]
    elif not variables['lipL']:
        # Rhombic dodecahedron (菱形十二面体) with square XY plane
        box = [[pbcx, 0, 0], [0, pbcy, 0], [0.5 * pbcx, 0.5 * pbcx, pbcz]]
    elif "hexagonal".startswith(options["-pbc"].value):
        box = [[pbcx, 0, 0], [math.sin(math.pi / 6) * pbcx, pbcy, 0], [0, 0, pbcz]]
    else:  # optimal packing; rhombic dodecahedron with hexagonal XY plane
        box = [[pbcx, 0, 0], [math.sin(math.pi / 6) * pbcx, pbcy, 0], [pbcx / 2, pbcy / 3, pbcz]]

    # Override lattice vectors if they were set explicitly
    box[0] = pbcSetX or box[0]
    box[1] = pbcSetY or box[1]
    box[2] = pbcSetZ or box[2]

    grobox = (box[0][0], box[1][1], box[2][2],
              box[0][1], box[0][2], box[1][0],
              box[1][2], box[2][0], box[2][1])

    pbcx, pbcy, pbcz = box[0][0], box[1][1], box[2][2]

    rx, ry, rz = pbcx + 1e-8, pbcy + 1e-8, pbcz + 1e-8

    pbc = pbcx,pbcy,pbcz
    rest_args['protein'] = protein
    rest_args['wholecenter'] = wholecenter
    rest_args['resi'] = resi
    rest_args['atid'] = atid
    rest_args['box'] = box
    rest_args['rx'] = rx
    rest_args['ry'] = ry
    rest_args['rz'] = rz
    rest_args['grobox'] = grobox


    return variables,options,lipids,pbc,rest_args