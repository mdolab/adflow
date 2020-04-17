
header_str = ''' pyf_processor.py is used to automatically process a
particular form of pyf file that allows the specification of both the
real and complex wrapping functions in one place. This is eliminates
duplication in the pyf files.

pyf_processor.py is called in the following manner:

python pyf_processor.py <real|complex> <pyf_file>

The resulting "clean" pyf file is written to pyf_file.autogen. This
.autogen file is what should be passed to f2py to generate the actual
wrapper.
'''
import sys

if len(sys.argv) != 3:
    print(header_str)
    sys.exit(1)
# end if

if sys.argv[1] in ['real','complex']:
    mode = sys.argv[1]
else:
    print(header_str)
    sys.exit(1)
# end if

try:
    f = open(sys.argv[2],'r')
    orig_lines = f.readlines()
except:
    print('Error opening/reading pyf file!')
    sys.exit(1)
# end try

# Open the new pyffile.
g = open(sys.argv[2] + '.autogen','w')

# Write the warning header
g.write('!'+'#'*78+'!'+'\n')
g.write('!' + ' '*27 + 'DO NOT MODIFY THIS FILE!'  + ' '*27 + '!' + '\n')
g.write('!' + ' '*27 + 'MODIFY adflow.pyf INSTEAD!'  + ' '*27 + '!' + '\n')
g.write('!'+'#'*78+'!'+'\n')

# Start going through the lines:
cur_mode = 'both'
for i in range(len(orig_lines)):
    if '#ifdef USE_COMPLEX' in orig_lines[i]:
        cur_mode = 'complex'
    elif '#ifndef USE_COMPLEX' in orig_lines[i]:
        cur_mode = 'real'
    elif '#else' in orig_lines[i]:
        # Flip the mode:
        if cur_mode == 'complex':
            cur_mode = 'real'
        elif cur_mode == 'real':
            cur_mode = 'complex'
        else:
            print('Error occured. Mismatched #else statement')
            sys.exit(1)
        # end if
    elif '#endif' in orig_lines[i]:
        cur_mode = 'both'
    else:
        # We have a normal line to write:
        if cur_mode == 'both' or cur_mode == mode:

            # We now have to check for real/integer types and process:
            if 'kind=inttype' in orig_lines[i]:
                orig_lines[i] = orig_lines[i].replace('kind=inttype','kind=4')
            if 'real(kind=realtype)' in orig_lines[i]:
                if mode == 'real':
                    orig_lines[i] = orig_lines[i].replace('kind=realtype','kind=8')
                else:
                    orig_lines[i] = orig_lines[i].replace(
                        'real(kind=realtype)','complex(kind=8)')
                # end if
            # end if
            g.write(orig_lines[i])
        # end if
    # end if
# end for

g.close()




