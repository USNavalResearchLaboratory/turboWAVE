# Program to translate OpenSCAD geometry files into turboWAVE equivalent.
# Not intended to be a general translator, supports limited subset of OpenSCAD.
# Contributed by G. Wijekoon.

print('Input File name. Ex: openscad or input help')
input2=input()
if input2=='help':
    he='''How to use:

	Enter input file name (not including extension). Input file must be txt.
	Program will then attempt to translate. It will print translation to console. The translation is after the line of ####.
	Enter y or n to accept translation. If y then it will create geometry.txt if it does not exist and write translation to it.
	Warning: It will overwrite everything in file.
	After writing to file, you will need to enter in the elements of the union named all.
	This is because Turbowave displays only one object, so all objects must be combined into one union.

	Other things you may have to do:

		If you have a difference, it will show up like below.
		You may or may not have to put the element name \'difference\' in union all

			new region intersection difference
			{
				elements = { mainbody , u0 }
			}

		If you have multiple differences, you will have to rename them because they all have the same name.

		If you have an object with a translation and rotate, you may have to change the order that they are shown in the geometry file.
		This is because TW rotates about global axis not object axis while openscad does the opposite.

	Input File Formatting:

		Object names:
		Must be put inline with Object. Must be last part in line. Format "//n [ObjectName]".
		No white space or percent signs are allowed in the name.

		Example OpenSCAD Code:
			union() { //n myunion
				cube(15, center=true); //n mycube
			}

	Brackets:

		No more than 1 bracket per line.
		Everything must be inside at least 1 set of brackets including top level functions.
		Example OpenSCAD Code:
			{
				sphere(10);
				cube(15.67, center=true);
			}
			// translate is not a child but still must be in brackets
			{
				translate([0,1,-0.5]) {
					cube([5,2,1]);
				}
			}

	Difference:

		When using a difference, all childs must be in union except for 1st
		Example OpenSCAD Code:
		{
			difference() {
				cube([20,10,5],center=true);  // <-- First child/object which is being cut
				// Objects being cut from first child
				union() {
					cylinder(2,r=3,center=true);
					translate([0,1,-0.5]) {
						cube([5,2,1]);
					}
					translate([-5,-3,-0.5]) {
						cube([5,2,1]);
					}
				}
			}
		}

	Indents and Spaces:

		Spaces are ok. Do not use notepad indents. Indents in openscad are fine.

	'''
    print(he)
    print('Input File name. Ex: openscad')
    input2=input()
    input1=input2+'.txt'
else:
    input1=input2+'.txt'

#Actual Code:

import re
def union(block,line,mods,nums,glob):
    '''This explains this function'''
    nl=line+1
    elements=[]
    l=1
    elements,glob,nums,nl = build(block,nl,nums,glob,l)
    elements= ' , '.join(elements)
    mods= '\n'.join(mods)
    un=nums[4]
    un=un+1
    nums[4]=un

    if (Section[block][line]).find('//n') != -1: # Detect if object was given name
        pos=(Section[block][line]).find('//n')+3 #find position of //n in line
        name=Section[block][line][pos:] #set name everything after //n
    else:
        name= 'u'+str(un-1)
    tp= 'new region union '+str(name)+'\n' '{' "\n" + '    elements = { '+str(elements) + " } \n    "+ str(mods) + "\n" '}'
    element=[]
    element.append(name)
    glob.append(tp)
    #nums= [0r,1s,2cy,3cn,4un,5inter,6diff]
    print('union break')
    return glob, nl , nums, element

def intersection(block,line,mods,nums,glob):
    nl=line+1
    elements=[]
    elements,glob,nums,nl = build(block,nl,nums,glob,1)
    print('int elements = ' + str(elements))
    elements= ' , '.join(elements)
    mods= '\n'.join(mods)
    inter=nums[5]
    if (Section[block][line]).find('//n') != -1: # Detect if object was given name
        pos=(Section[block][line]).find('//n')+3 #find position of //n in line
        name=Section[block][line][pos:] #set name everything after //n
    else:
        name= 'inter'+str(inter)
    tp= 'new region intersection '+ str(name) +'\n' '{' "\n" + '    elements = { '+str(elements) + " } \n    "+ str((mods)) + "\n" '}'
    inter=inter+1
    element=[]
    element.append(name)
    glob.append(tp)
    nums[5]=inter
    #nums= [0r,1s,2cy,3cn,4un,5inter,6diff]
    return glob, nl , nums , element


## Modifier Functions
def translate(block,nl):
    TC=  [float(d) for d in re.findall(r'-?\d+\.?\d*',Section[block][nl])]
    dx=str(TC[0])
    dy=str(TC[1])
    dz=str(TC[2])
    Tm= '    translation = ( ' + dx +' , ' +dy+ ' , ' +dz + ') '
    return Tm

def rotate(block,nl):
    RC= [float(d) for d in re.findall(r'-?\d+\.?\d*',Section[block][nl])]
    rx=str(RC[0])
    ry=str(RC[1])
    rz=str(RC[2])
    Rx = '    rotate about x =  ' + rx + '\n'
    Ry = '    rotate about y  =  ' + ry + '\n'
    Rz = '    rotate about z  =  ' + rz
    if RC[0]==0:
        Rx = ''
    if RC[1]==0:
        Ry = ''
    if RC[2]==0:
        Rz = ''
    Rm= Rx + Ry + Rz
    return Rm
    return Rm
######################################
def Cube(block,line,mods,rec):
    #Cube_coords=[]
    CoordsSize= re.findall(r'\d+\.?\d*',Section[block][line]) #get values in the line
    cs= 0 # clear cs variable every loop
    if len(CoordsSize) == 1: #if actual cube ie 1 number given
        cs= float(CoordsSize[0])
        if (Section[block][line]).find('center=true')!=-1:  # check if center is true
            x0=  -cs/2
            x1= cs/2
            y0= -cs/2
            y1= cs/2
            z0= -cs/2
            z1= cs/2
        else:
            x0= 0
            x1= cs
            y0= 0
            y1= cs
            z0= 0
            z1= cs
    else: # if not cube but rect with xyz dimensions
        cs= CoordsSize
        cs=[float(i) for i in cs]
        print(cs)
        if (Section[block][line]).find('center=true')!=-1:  # check if center is true
            x0= -cs[0]/2
            x1= cs[0]/2
            y0= -cs[1]/2
            y1= cs[1]/2
            z0= -cs[2]/2
            z1= cs[2]/2
        else:
            x0= 0
            x1= cs[0]
            y0= 0
            y1= cs[1]
            z0= 0
            z1= cs[2]
    if (Section[block][line]).find('//n') != -1: # Detect if object was given name
        pos=(Section[block][line]).find('//n')+3 #find position of //n in line
        name=Section[block][line][pos:] #set name everything after //n
    else:
        rec=nums[0]
        name= 'r'+str(rec)


    modifications='\n'.join(mods)
    print('cube found name: '+name+' block:'+ str(block) +' line:'+str(line))
    Cube_coords= ' ( ' + str(x0)+ ' , '+ str(x1) + ' , ' +  str(y0)+ ' , '+ str(y1) + ' , ' + str(z0)+ ' , '+ str(z1) + ' ) '
    Cube= 'new region rect ' + str(name) + '\n' '{ \n' '    ' 'bounds = '+ Cube_coords  + '\n' + modifications+ '\n }'
    return Cube,name

###############################################3
def build(block,nl,nums,glob,l):
    elements=[]
    mods=[]

    print('new build started', mods)
    #l=1
    r=0
    diff=0
    tro=-1
    tror=0
    trig=0
    ml=-1 #mod list

    while True:# l != r and l != 0 :
            print('L'+  str(nl))
            #print(nl)
            if nl >= len(Section[block])-1:# or (Section[block][nl]).find('}') == 0:
                #print(len(Section[block]))
                print('len break')
                break
            r=(Section[block][nl]).count('}')+r
            l=(Section[block][nl]).count('{')+l
            print('r='+str(r))
            print('l='+str(l))
            if l == r and l!=0:
                print('break')
                break
            #
            # Modifiers /non-objects
            #modifiers translate and rotate. When found, adds modifier to list: mods
            if  (Section[block][nl]).find('translate') == 0:
                print('translate found block:'+ str(block) +' line:'+str(nl))
                e=translate(block,nl)
                ml=ml+1
                mods.append(e)
                tro=0
            if  (Section[block][nl]).find('rotate') == 0:
                print('rotate found block:'+ str(block) +' line:'+str(nl))
                e=rotate(block,nl)
                mods.append(e)
                tro=0
                ml=ml+1
            # used to remove modifier from mods when loop is out of that modifier's nest
            if tro != -1:

                tro=tro + (Section[block][nl]).count('{')
                tror=(Section[block][nl]).count('}')+tror
                print(tro, tror)
                print(ml)
                if tro== tror:
                    del mods[len(mods)-1]
                    ml=ml-1
                    if ml != -1:
                        print('chain not done')
                        tro=tro+1
                    else:
                        tro=-1
                        tror=0
                        print('chain done')
                    print(mods)
            ################################## End modifiers

            if (Section[block][nl]).find('difference') == 0:
                print('difference found')
                di= "complement = true   "
                diff=2
                trig=1


            # Objects with childs


            if (Section[block][nl]).find('union') == 0:
                print('union found block:'+ str(block) +' line:'+str(nl))
                glob, nl, nums,element = union(block,nl,mods,nums,glob)
                element= ' , '.join(element)
                elements.append(element)
                if diff==2:
                    diff=1

            if (Section[block][nl]).find('intersection') == 0:
                print('intersection found block:'+ str(block) +' line:'+str(nl))
                glob, nl, nums,element = intersection(block,nl,mods,nums,glob)
                element= ' , '.join(element)
                elements.append(element)
                if diff==2:
                    diff=1

            # Objects without childs
            if (Section[block][nl]).find('sphere') == 0:
                        s=nums[1]
                        s=s+1
                        nums[1]=s
                        if (Section[block][nl]).find('//n') != -1: # Detect if object was given name
                            pos=(Section[block][nl]).find('//n')+3 #find position of //n in line
                            name=Section[block][nl][pos:] #set name everything after //n
                        else:
                            name= 's'+str(s-1)
                        elements.append(name)
                        Sphere_radius = re.findall(r'\d+', Section[block][nl])
                        sr = float(Sphere_radius[0]) #turn radius value in list into float

                        Sphere= 'new region true_sphere ' + str(name) + '\n' '{' "\n" '    ' 'radius = '+ str(sr) + "\n" + str(mods)[2:-2] + '\n }'
                        glob.append(Sphere)
                        if diff==2:
                            diff=1

                        print('sphere found block:'+ str(block) +' line:'+str(nl))
                        #nums= [0re,1s,2cy,3cn,4un,5inter,6diff]

            if (Section[block][nl]).find('cylinder') == 0:
                if diff==2:
                    diff=1
                cc= re.findall(r'\d+',Section[block][nl])
                temp= cc.copy() # make temp list for checking if r1=r2
                del temp[0] # remove height in temp
                temp=set(temp)# remove duplicate values, if r1=r2 then list size = 1
                modifications='\n'.join(mods)
                if len(cc)==2 or len(temp)==1: # check if only height and radius given or if temp list size=1
                    cy=nums[2]
                    if (Section[block][nl]).find('//n') != -1: # Detect if object was given name
                        pos=(Section[block][nl]).find('//n')+3 #find position of //n in line
                        name=Section[block][nl][pos:] #set name everything after //n
                    else:
                        name= 'cy'+str(cy)
                    print('cylinder found'+str(name))
                    elements.append(name)
                    cy=cy+1
                    nums[2]=cy
                    Cylinder = 'new region cylinder '+str(name) + '\n' '{' "\n" '    ' 'length = '+ str(cc[0]) + "\n" '    ' 'radius = '+ str(cc[1])+ "\n" +modifications +'\n }'
                    glob.append(Cylinder)
                else:
                    cn= nums[3]
                    if (Section[block][nl]).find('//n') != -1: # Detect if object was given name
                        pos=(Section[block][nl]).find('//n')+3 #find position of //n in line
                        name=Section[block][nl][pos:] #set name everything after //n
                    else:
                        name= 'cn'+str(cn)
                    if (Section[block][nl]).find('center=true')!=-1:
                        delta= 'translation = ( 0 , 0 ,'+ str(int(cc[0])/2) +')'
                        mods.append(delta)
                    elements.append(name)
                    Cone = 'new region cone '+str(name) + '\n' '{' "\n" '    ' 'length = '+ str(cc[0]) + "\n" '    ' 'tip radius = '+ str(cc[2])+ "\n" '    ' 'base radius = '+str(cc[1])+ "\n" +modifications+  '\n}'
                    cn=cn+1
                    nums[3]=cn
                    del mods[(len(mods)-1)]
                    glob.append(Cone)



            if (Section[block][nl]).find('cube') == 0:
                rec=nums[0]
                tp,name=Cube(block,nl,mods,rec)
                elements.append(name)
                rec=rec+1
                nums[0]=rec
                if diff==2:
                    diff=1
                glob.append(tp)
                print(elements)

            nl=nl+1
            if diff==1:
                mods.append(di)
                diff=0
            if nl >= len(Section[block])-1:# or (Section[block][nl]).find('}') == 0:
                #print(len(Section[block]))
                print('len break')
                break

    if trig==1:
        elements= ' , '.join(elements)
        tp= 'new region intersection difference \n' '{' "\n" + '    elements = { '+str(elements) + " } \n} "
        glob.append(tp)
    return elements,glob,nums,nl


Input = open(input1, 'r')
Output = open('geometry.txt','w')
std= Output
scad = Input.read()


########################
#Pt1: Detect big blocks
#Function code
#notes: all functions must come first. Any seperate objects must be last.

text = scad
#Change text formatting
for line in text: # remove all uppercase
       line1 = text.lower()
line3 = [line.replace('\t' and '\t\t' and '	', ' ') for line in line1] # remove indents
lines = ''.join([line.replace(' ', '') for line in line3]) # remove spaces



#Detects big nests
blocks=[]
snum=0 #section number
r=0 # number of } brackets
l=0 # number of { brackets
lin=0
temp=0

#What this does: goes through every line and counts { (variable l) and } (variable r)
#                brackets, when l = r, creates list with all lines from first block
#                as the first entry. Then removes the block from lines (list).
#print(lines)
for space in lines:
    r=space.count('}')+r
    l=space.count('{')+l
    #print(l,r)
    if l == r and temp != r: #section found
    #    print('section found')
        temp=r #prevent making new section when l=r and there are blank lines
        space= lines[:lin] # defines section
        c= space.count('') #counts number of characters
        se= lines[:c] #creates string with only the first section
        lines=lines[c:] #redefines lines by excluding the first section
        blocks.append(se) # adds first section to Sections list
        snum=snum+1 #count number of sections for later
    lin=lin+1
   # now when it loops, the 2nd section is now the first section
del temp
Section=[]
# Splits each section into multiple lines. Then makes each line a part of the list (2d list)
for block in range(0,snum,1): #for each section
    part=(((blocks[block]).splitlines())) # split each line of each section
    Section.append(part)
del block
del blocks
#pprint(Section)

###################################
#Pt2: Builds list (glob) with all text to write to file
#Function Code

rec=0 #cube number
s=0 #sphere number
cy=0 #cylinder number
cn=0 #cone number
un=0 #union number
inter=0 #intersection number
diff=0 #difference number
nums= [rec,s,cy,cn,un,inter,diff]
elements=[]
l=-1
r=0
glob=[] # "global" list that carries everything to be written to file
nl=0
block=0
line=0
g=[]

for block in range(0,len(Section),1):
    elements,glob,nums,nl=build(block,0,nums,glob,l)
    g.extend(glob)
    g.append('//Section End')
    glob=[]
    mods=[]
    elements=[]
    l=-1
    r=0


all= '// Modify all here  \n////////////////////////////////////////// \n' + str('new region union all \n{\n    elements = {  } \n}')
g.append(all)

glb='\n'.join(g)
print('################################################')
#glb= [line.replace("['", '') for line in glb]
print(glb)

print('Write to file? y or n')
if input()=='y':
   #Write geometry to file
    std.write(glb)
    print('Geometry written')
else:
   print('Geometry not written')
std.close()

#######################################
