
import sys
import pandas as pd
import numpy as np

# set this to change the precision of the values printed, it must be a string
DIGITS = '6'

def write_bngl_function(input_fname, param):
    # load the input file, make sure that the values are floats
    df = pd.read_csv(
        input_fname, delim_whitespace=True, comment='#', 
        names=['time', 'value'], dtype={'time':np.float64, 'value':np.float64})
                                        
    #print(df['time'])
    #print(df['value'])
    
    if df['time'][0] != 0:
        print("First time must be 0")
        sys.exit(1)
    
    ind = "    "
    for i in range(1, len(df)):
        print(
            ind +
            "(if((" + param + "<" + str.format('{0:.' + DIGITS + '}', df['time'][i]) + "), " + 
            "(" + str.format('{0:.' + DIGITS + '}', df['value'][i - 1]) + "), \\") 

    # last value
    print(ind + "  (" + str.format('{0:.' + DIGITS + '}', df['value'].iloc[-1]) + ") \\") 

    print(ind, end='') 
    for i in range(1, len(df)):
        print("))", end = '')
    print()

def main():
    if len(sys.argv) != 3:
        print("Expecting input file as 1st argument and input variable name (usually time) as 2nd argument")
        sys.exit(1)
    
    write_bngl_function(sys.argv[1], sys.argv[2])
    

if __name__ == "__main__":
    main()
    