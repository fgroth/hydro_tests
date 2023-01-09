#!/usr/bin/python
# author: Frederick Groth (fgroth@usm.lmu.de)
#
# aim:
#  can be used to standardize parameter files or convert them between different codes.
#  compare also param_conversion.param, contains header (%%%), followed by parameters
#  in desired order including comments
#  header format: '<code_name>:<#symbols>' in the same order as columns
#
# usage:
#  python standardize_param.py <input parameter file>
# in this case, both as input and output type OpenGadget3 is assumed
# or 
#  python standardize_param.py <input parameter file> <input code type> <output code type>
# or 
#  python standardize_param.py <input parameter file> <input code type> <output code type> <mode>
# where mode is a bitwise read int:
#  1 : ask for deletion (0 is not asking)
#  2 : append or delete by default (0 is appending)
#  4 : ask for values directly (0 is not asking)
# a file with the name <input parameter file>_new is produced, having a standardized format.
# if necessary, missing values have to be added.

import sys

if len(sys.argv) > 1:
    file = sys.argv[1]
else:
    print('filename as input argument missing')
    exit()

if len(sys.argv) > 3:
    input_code_type = sys.argv[2]
    output_code_type = sys.argv[3]
else:
    input_code_type="OpenGadget3"
    output_code_type="OpenGadget3"

if len(sys.argv) > 4:
    mode = int(sys.argv[4])
else:
    mode=0

columns = 38

template='param_conversion.param'
new_file = file+'_new'
f_new = open(new_file, 'w')

comment_chars = ["%", ";", "#", "!"]

def read_parameters(file):
    """
    read in parameters from file into dictionary

    input parameter:
      file: String, filename
    
    return value:
      parameters: dictionary of the form {"parameter_name":"value"}
    """
    parameters={}
    # read original parameter file
    with open(file) as f:
        for line in f:
            for char in comment_chars:
                line = line.split(char)[0] #remove comments
            line=line.split('\n')[0] #remove linebreaks
            if len(line.split()) != 0: # line not empty (ignoring spaces)
                parameters[line.split()[0]]=line.split()[-1] # add the parameter to the Dictionary parameters
    f.close()
    return parameters

def read_header(line,codes,positions):
    """
    read one line the header of template file (param_conversion.py) and add the information to the output lists.

    input parameters:
      line : String containing a header line
     output parameters:
       codes: list of Strings (names of the different codes for which parameters are present)
       positions: list of Integers (length)
    
    return value:
      0
    """
    new_line = line[4:]
    new_line=new_line.split('\n')[0] #remove linebreaks
    info = new_line.split(":")
    code = info[0]
    nchar = int(info[1])
    codes.append(code)
    positions += [positions[-1]+nchar]
    return 0

def get_parameters(line,code,codes,positions):
    """
    get name of parameter

    input parameters:
      line: String, line from template file (param_conversion.param)
      code: String, code for which the parameter name should be found.
      codes: list of Strings, codes as defined in header of template file.
      positions: list of Integers, length of blocks within template file (compare header)
    
    return value:
      param: list of String(s), name(s) of parameter within desired code.
    """
    i_code = codes.index(code)
    param = line[positions[i_code]:positions[i_code+1]]
    # remove spaces
    param = param.replace(" ","")
    return param # return name and alternative name in one string, as in template file.
    
def add_line(line,new_file,parameters={}):
    """
    write a line of the parameter into the new, standardized parameter file.

    input parameters:
      line: String, line within original parameter file
      new_file: String, filename of new file that the line should be added to.
    
    return value:
      0
    """
    for char in ['%']: #comment_chars:
        new_line=line.split(char)[0] #remove comments
    new_line=new_line.split('\n')[0] #remove linebreaks
    if len(new_line) != 0:
        value=''
        # only take the columns that we need
        input_param = get_parameters(new_line,input_code_type,codes,positions) # this is String and might contain an alternative names separated by "/".
        output_param = get_parameters(new_line,output_code_type,codes,positions).split("/")[0] # always use first (prefered) name
        # check if parameter exists in the desired code version
        if len(input_param) != 0 and len(output_param) == 0:
            if mode >> 1 & 1:
                tmp = "n"
            else:
                tmp = "y"
            print(input_param+" exists only in "+input_code_type)
            if mode >> 0 & 1:
                print("Delete it (y) ar append it to the end (n)? Default is "+tmp)
                tmp2 = ""
                while (tmp2 not in ["y","n"]):
                    tmp2 = raw_input("enter y/n : ")
                    if len(tmp2) == 0:
                        tmp2 = tmp
                    if tmp2 == "y":
                        for param in input_param.split("/"):
                            if param in parameters:
                                del parameters[param]
            else:
                if tmp == "y":
                    for param in input_param.split("/"):
                            if param in parameters:
                                del parameters[param]
                    print("Parameter was deleted")
                else:
                    print("Parameter is appened")
            return 1 # not appended yet
        elif len(input_param) == 0 and len(output_param) != 0:
            print("WARNING: "+output_param+" exists only in "+output_code_type+" not in "+input_code_type)
            if mode >> 2 & 1:
                value = raw_input("Value: ")
        elif len(input_param) != 0 and len(output_param) != 0: # parameter exists
            # now we can assign values
            value_present = False
            for param in input_param.split("/"):
                if param in parameters:
                    if not value_present: # only read first occurence
                        value = parameters[param]
                        del parameters[param]
                        value_present = True
            if not value_present:
                print('WARNING: No value for '+input_param+' provided')
        while len(output_param) < columns:
            output_param+=' ' #make it nice looking
        output_line=output_param+value+'\n'
        f_new.write(output_line)
    else:
        f_new.write(line)
    return 0

parameters = read_parameters(file)


#add parameters
codes = []
positions = [0]
with open(template) as f:
    for line in f:
        # read header
        if line.startswith("%%%"):
            read_header(line, codes,positions)
        else:
            add_line(line,f_new,parameters=parameters)
f.close()
            
#go over elements left in list
for line in parameters:
    value = parameters[line]
    while len(line) < columns:
        line+=' ' #make it nice looking
    line+=value
    line+='\n'        
    f_new.write(line)


f_new.close()
