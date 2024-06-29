import sys

def output(variables,options,lipids,pbc,rest_args):
    #Write the output
    slen = rest_args['solvent'] and len(rest_args['sol']) or 0
    print("; NDX Solvent %d %d" % (1+rest_args['plen']+rest_args['mlen'], rest_args and rest_args['plen']+rest_args['mlen']+slen or 0))
    print("; NDX System %d %d" % (1, rest_args['plen']+rest_args['mlen']+slen))

    # Open the output stream
    oStream = options["-o"] and open(options["-o"].value,"w") or sys.stdout

    # Print the title
    if rest_args['membrane'].atoms:
        title  = "INSANE! Membrane UpperLeaflet>"+":".join(variables['lipU'])+"="+":".join([str(i) for i in variables['numU']])
        # title += " LowerLeaflet>"+":".join(lipL)+"="+":".join([str(i) for i in numL])

        if rest_args['protein']:
            title = "Protein in " + title
    else:
        title = "Insanely solvated protein."

    print(title,file = oStream)

    # Print the number of atoms
    print("%5d"%(len(rest_args['membrane'])+len(rest_args['sol'])),file = oStream)

    # Print the atoms
    id = 1
    if rest_args['protein']:
        for i in range(len(rest_args['protein'])):
            at,rn,ri = rest_args['protein'].atoms[i][:3]
            x,y,z    = rest_args['protein'].coord[i]
            oStream.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(ri,rn,at,id,x,y,z))
            id += 1
    if rest_args['membrane']:
        for i in range(len(rest_args['membrane'])):
            at,rn,ri = rest_args['membrane'].atoms[i][:3]
            x,y,z    = rest_args['membrane'].coord[i]
            oStream.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(ri,rn,at,id,x,y,z))
            id += 1
    if rest_args['sol']:
        # Print the solvent
        print("\n".join([i[0]+"%8.3f%8.3f%8.3f"%i[1] for i in rest_args['sol']]),file = oStream)

    print("%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n"%rest_args['grobox'],file = oStream)


    if options["-p"].value:
         # Write a rudimentary topology file
         top = open(options["-p"].value,"w")
         print('#include "martini.itp"\n',file=top)
         print('[ system ]\n; name\n%s\n\n[ molecules ]\n; name  number'%title,file=top)
         if rest_args['']:
             print("%-10s %5d"%("Protein",1),top)
         print("\n".join("%-10s %7d"%i for i in rest_args['molecules']),file=top)
         top.close()
    else:
         print("\n".join("%-10s %7d"%i for i in rest_args['molecules']),file=sys.stderr)

    print("; This script, written by by Yuan He and Yuxuan Wang in Dr. Xubo Lin's group@Beihang University, can be used to generate lipid nanobubble.")
