# Because for some god forsaken reason julia doesn't have a switch-case struct
# by default we use a package to give it to us. It's overkill but whatever.
using Match

# The parsing of the state function is significantly easier if we run the julia
# code and thus define the state structures into memory.
# Honestly this more or less boils down to getting number of dimensions.
# This also ensures that there are no obvious bugs in your julia code.
# include("structures_test.jl")

# Actually, scratch that. I'll just make it an input argument for the parsing.

# first we generate the header file.
function generate_header(file)
    ioc = open("pfc_header.h","w+")    # generate header file
    write(ioc, "// dependencies\n#include <mpi.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <hdf5.h>
#include <gsl/gsl_rng.h>
#define PI 2.0*acos(0.0)\n")    # write standard set of dependencies
    write_state_header(file, ioc)
    write_function_header(file, ioc)
    # write_io_header(file, ioc)  # this one will be a pain in the ass...
    close(ioc)
end

# parse and write the state in the header
function write_state_header(file, ioc)
    write(ioc, "typedef struct\n{\n")   # open the state structure
    ioj = open(file, "r")   # open julia file in write only mode
    while(!eof(ioj))
        lines = readuntil(ioj, "end\n")
        if startswith(lstrip(lines), "mutable struct")
            empty, lines = split(lines, "mutable struct ")  # strip away the struct
            type, lines... = split(chomp(lines),"\n")   # seperate out the type definition
            # match the type definition to appropriate c types
            clines = @match type begin
                "ReArrays"      =>  replace.(lstrip.(lines), r"(?<var>\w+)" => s"double* \g<var>;\n")
                "CmxArrays"     =>  replace.(lstrip.(lines), r"(?<var>\w+)" => s"fftw_complex* \g<var>;\n")
                "halfReArrays"  =>  replace.(lstrip.(lines), r"(?<var>\w+)" => s"double* \g<var>;\n")
                "dims"          =>  replace.(lstrip.(lines), r"(?<var>\w+)" => s"ptrdiff_t \g<var>;\n")
                "IntParams"     =>  replace.(lstrip.(lines), r"(?<var>\w+)" => s"int \g<var>;\n")
                "doubleParams"  =>  replace.(lstrip.(lines), r"(?<var>\w+)" => s"double \g<var>;\n")
                "transforms"    =>  replace.(lstrip.(lines), r"(?<var>\w+)" => s"fftw_plan \g<var>;\n")
                _ => ""
            end
            write(ioc, join(clines)*"\n")
        end
    end
    close(ioj)
    # add the parts of the structure that are purely for mpi codes
    write(ioc, lstrip("ptrdiff_t local_n0;       // local endpoint of a processor
ptrdiff_t local_0_start;  // real start of a processor
ptrdiff_t local_n1;       // Used for transposed fftw
ptrdiff_t local_1_start;  // Same.\n"))
    # add the rng
    write(ioc, "gsl_rng* rng;\n")
    # close the state structure
    write(ioc, "} state;\n\n")
end

# parse and write functions in the header
function write_function_header(file, ioc)
    # here we assume that non typecast functions and variables are of type void
    ioj = open(file, "r")   # open julia file in write only mode
    while(!eof(ioj))
        lines = readuntil(ioj, "end #function\n")
        # if startswith(lstrip(lines), "function")
        empty, lines = split(lines, "function ")  # strip away the function definition
        name, body... = split(chomp(lines),"\n")   # seperate out the function name from the body
        # here we only  care about the name since its just the header
        name, functype = split(name,")")
        name, args = split(name,"(")
        # deal with typecasting of the function
        name = @match functype begin
            "state"     =>  replace(name, r"(?<funcname>\w+)" => s"state* \g<funcname>(")
            "Int"       =>  replace(name, r"(?<funcname>\w+)" => s"int \g<funcname>(")
            "Float64"   =>  replace(name, r"(?<funcname>\w+)" => s"double \g<funcname>(")
            _           =>  replace(name, r"(?<funcname>\w+)" => s"void \g<funcname>(")
        end
        # split up any arguments
        args = split(args, ",")
        # deal with typecasting
        for i in eachindex(args)
            arg, type = split(args[i], "::")
            args[i] = @match type begin
                "state"     =>  replace(arg, r"(?<var>\w+)" => s"state* \g<var>,")
                "Int"       =>  replace(arg, r"(?<var>\w+)" => s"int \g<var>,")
                "Float64"   =>  replace(arg, r"(?<var>\w+)" => s"double \g<var>,")
                _           =>  replace(arg, r"(?<var>\w+)" => s"void \g<var>,")
            end
            if i == length(args)
                args[i] = chop(args[i])*");\n"
            end
        end
        args = join(args)
        cline = name*args
        write(ioc, cline)
        # end
    end
    close(ioj)
end

# parse and write functions in the header
function write_functions(file, dims, ioc)
    # here we assume that non typecast functions and variables are of type void
    ioj = open(file, "r")   # open julia file in write only mode
    while(!eof(ioj))
        lines = readuntil(ioj, "end #function\n")
        prefunc, lines = split(lines, "function")
        # "function"*lines
        # if startswith(lstrip(lines), "function")
            # empty, lines = split(lines, "function ")  # strip away the function definition
            name, body... = split(chomp(lines),"\n",keepempty=false)   # seperate out the function name from the body
            # println(body)
            bodycopy = body
            # here we only  care about the name since its just the header
            name, functype = split(name,")")
            name, args = split(name,"(")
            # deal with typecasting of the function
            name = @match functype begin
                "state"     =>  replace(name, r"(?<funcname>(\w|\d)+)" => s"state* \g<funcname>(")
                "Int"       =>  replace(name, r"(?<funcname>(\w|\d)+)" => s"int \g<funcname>(")
                "Float64"   =>  replace(name, r"(?<funcname>(\w|\d)+)" => s"double \g<funcname>(")
                _           =>  replace(name, r"(?<funcname>(\w|\d)+)" => s"void \g<funcname>(")
            end
            # split up any arguments
            args = split(args, ",")
            cargs = similar(args)
            # deal with typecasting
            for i in eachindex(args)
                arg, type = split(args[i], "::")
                cargs[i] = @match type begin
                    "state"     =>  replace(arg, r"(?<var>(\w|\d)+)" => s"state* \g<var>,")
                    "Int"       =>  replace(arg, r"(?<var>(\w|\d)+)" => s"int \g<var>,")
                    "Float64"   =>  replace(arg, r"(?<var>(\w|\d)+)" => s"double \g<var>,")
                    _           =>  replace(arg, r"(?<var>(\w|\d)+)" => s"void \g<var>,")
                end
                args[i] = @match type begin
                    "state"     =>  replace(arg, r"(?<var>(\w|\d)+)" => s"\g<var>")
                    "Int"       =>  replace(arg, r"(?<var>(\w|\d)+)" => s"\g<var>")
                    "Float64"   =>  replace(arg, r"(?<var>(\w|\d)+)" => s"\g<var>")
                    _           =>  replace(arg, r"(?<var>(\w|\d)+)" => s"\g<var>")
                end
                if i == length(args)
                    cargs[i] = chop(cargs[i])*")\n"
                end
            end

            args2 = join(cargs)
            cline = name*args2*"{\n"
            write(ioc, cline)


            parse_local_function_variables(bodycopy, args, ioc, 1)
            # println(body)
            # println(bodycopy.*"\n")
            # println(body.*"\n")
            parse_function_body(body, args, dims, ioc)
        # end
    end
    close(ioj)
end

function generate_state(file, D)
    ioc = open("state.c","w+")    # generate state file
    write(ioc, "// Define dependencies
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include \"pfc_header.h\"\n")    # write standard set of dependencies
    # write_state_functions(file, ioc)

    # the state file will always have 2 base functions
    # 1) creation
    # 2) destruction

    # creation
    @match D begin
        2 => write(ioc, "state* create_state(Nx, Ny)\n{\n")
        3 => write(ioc, "state* create_state(Nx, Ny, Nz)\n{\n")
        _ => error("Dimensionality is assumed to be 2 or 3. Check your code.")
    end
    write(ioc, "ptrdiff_t local_alloc;\n")
    write(ioc, "state* s = malloc(sizeof(state));\n")
    write(ioc, "if (s == NULL) return NULL;\n")
    # Allocate memory based on number of dimensions
    @match D begin
        2 => write(ioc, "local_alloc = fftw_mpi_local_size_2d_transposed(Nx,(Ny>>1)+1, MPI_COMM_WORLD,
                                        &s->local_n0, &s->local_0_start,
                                        &s->local_n1, &s->local_1_start);\n")
        3 => write(ioc, "local_alloc = fftw_mpi_local_size_3d_transposed(Nx,Ny,(Nz>>1)+1, MPI_COMM_WORLD,
                                        &s->local_n0, &s->local_0_start,
                                        &s->local_n1, &s->local_1_start);\n")
        _ => error("Dimensionality is assumed to be 2 or 3. Check your code.")
    end
    # iterate through the julia file and allocate the state arrays
    ioj = open(file, "r")   # open julia file in read only mode
    while(!eof(ioj))
        lines = readuntil(ioj, "end\n")
        if startswith(lstrip(lines), "mutable struct")
            empty, lines = split(lines, "mutable struct ")  # strip away the struct
            type, lines... = split(chomp(lines),"\n")   # seperate out the type definition
            # match the type definition to appropriate c types
            clines = @match type begin
                "ReArrays"      =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"s->\g<var> = fftw_alloc_real(2*local_alloc);\n")
                "CmxArrays"     =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"s->\g<var> = fftw_alloc_complex(local_alloc);\n")
                "halfReArrays"  =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"s->\g<var> = fftw_alloc_real(local_alloc);\n")
                "dims"          =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"s->\g<var> = \g<var>;\n")
                "IntParams"     =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"s->\g<var> = 0;\n")
                "doubleParams"  =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"s->\g<var> = 0.0;\n")
                _ => ""
            end
            write(ioc, join(clines)*"\n")
        end
    end
    # error checking to ensure arrays have been allocated
    seekstart(ioj)
    write(ioc, "if (")
    while(!eof(ioj))
        lines = readuntil(ioj, "end\n")
        if startswith(lstrip(lines), "mutable struct")
            empty, lines = split(lines, "mutable struct ")  # strip away the struct
            type, lines... = split(chomp(lines),"\n")   # seperate out the type definition
            # match the type definition to appropriate c types
            clines = @match type begin
                "ReArrays"      =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"s->\g<var> == NULL ||\n")
                "CmxArrays"     =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"s->\g<var> == NULL ||\n")
                "halfReArrays"  =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"s->\g<var> == NULL ||\n")
                _ => ""
            end
            write(ioc, join(clines)*"\n")
        end
    end
    write(ioc, "false)\n{\n")
    seekstart(ioj)
    while(!eof(ioj))
        lines = readuntil(ioj, "end\n")
        if startswith(lstrip(lines), "mutable struct")
            empty, lines = split(lines, "mutable struct ")  # strip away the struct
            type, lines... = split(chomp(lines),"\n")   # seperate out the type definition
            # match the type definition to appropriate c types
            clines = @match type begin
                "ReArrays"      =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"free(s->\g<var>);\n")
                "CmxArrays"     =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"free(s->\g<var>);\n")
                "halfReArrays"  =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"free(s->\g<var>);\n")
                _ => ""
            end
            write(ioc, join(clines)*"\n")
        end
    end
    write(ioc, "return NULL;\n}\n")
    seekstart(ioj)
    # write the fftw plans. This is somewhat more complicated...
    # going to just assume that we have one forward and one back transform.
    # order will be assumed to be forward first
    re4plan = nothing
    im4plan = nothing
    fname = nothing
    bname = nothing
    while(!eof(ioj))
        lines = readuntil(ioj, "end\n")
        if startswith(lstrip(lines), "mutable struct")
            empty, lines = split(lines, "mutable struct ")  # strip away the struct
            type, lines... = split(chomp(lines),"\n")   # seperate out the type definition
            # match the type definition to appropriate c types
            if (type == "ReArrays" && re4plan == nothing)
                re4plan = lstrip(lines[1])
            end
            if (type == "CmxArrays" && im4plan == nothing)
                im4plan = lstrip(lines[1])
            end
            if (type == "transforms" && fname == nothing && bname == nothing)
                fname = lstrip(lines[1])
                bname = lstrip(lines[2])
            end
        end
    end
    @match D begin
        2 => begin
                write(ioc, "s->$fname = fftw_mpi_plan_dft_r2c_2d(Nx, Ny, s->$re4plan, s->$im4plan, MPI_COMM_WORLD,
                                                  FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);\n")
                write(ioc, "s->$bname = fftw_mpi_plan_dft_c2r_2d(Nx, Ny, s->$im4plan, s->$re4plan, MPI_COMM_WORLD,
                                                  FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);\n")
              end
        3 => begin
                write(ioc, "s->$fname = fftw_mpi_plan_dft_r2c_3d(Nx, Ny, Nz, s->$re4plan, s->$im4plan, MPI_COMM_WORLD,
                                                  FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);\n")
                write(ioc, "s->$bname = fftw_mpi_plan_dft_c2r_3d(Nx, Ny, Nz, s->$im4plan, s->$re4plan, MPI_COMM_WORLD,
                                                  FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);\n")
              end
        _ => error("Dimensionality is assumed to be 2 or 3. Check your code.")
    end
    # write the rng
    seekstart(ioj)
    write(ioc, "// Nathan style entropy
s->rng = gsl_rng_alloc (gsl_rng_default);
if (s->rng == NULL)
{
gsl_rng_free (s->rng);
return NULL;
}
int rank;
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
gsl_rng_set (s->rng, (
    (unsigned long)time(NULL)
     + (unsigned long)clock()
     + (unsigned long)getpid()
     + (unsigned long)getppid()) * (rank + 1));

// Barrier to avoid pipeline issues.
MPI_Barrier(MPI_COMM_WORLD);

return s;\n}\n")

    # destruction
    write(ioc, "void destroy_state(state* s)
{
  // Free memory
  if (s != NULL)
  {
  fftw_destroy_plan (s->$fname);
  fftw_destroy_plan (s->$bname);\n")
    while(!eof(ioj))
        lines = readuntil(ioj, "end\n")
        if startswith(lstrip(lines), "mutable struct")
          empty, lines = split(lines, "mutable struct ")  # strip away the struct
          type, lines... = split(chomp(lines),"\n")   # seperate out the type definition
          # match the type definition to appropriate c types
          clines = @match type begin
              "ReArrays"      =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"free(s->\g<var>);\n")
              "CmxArrays"     =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"free(s->\g<var>);\n")
              "halfReArrays"  =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"free(s->\g<var>);\n")
              _ => ""
          end
          write(ioc, join(clines)*"\n")
        end
    end
    write(ioc,"free(s);\n}\n\n}")

    close(ioj)
    close(ioc)

    return nothing
end

function generate_io(file)
    # io will generate a verbose save function that outputs all possible real fields.
    # it will then be up to the user to comment out the fields they don't want to save.
    ioc = open("io.c","w+")
    # define requirements
    write(ioc, "#include <mpi.h>        // MPI runtime
#include <hdf5.h>       // HDF5 header
#include <stdlib.h>     // malloc etc
#include <unistd.h>     // access() to check if file exists
#include <stdbool.h>    // Booleans
#include <assert.h>     // Assertions
#include <string.h>     // String manipulation sprintf etc...
#include <math.h>       // Needed for sqrt

#include \"pfc_header.h\"    // local libbinary header
#include \"io_convenience.c\"

#define FAIL -1
#define LEN 8

const int io_verbose = true;\n")

    # write the save function
    ioj = open(file, "r")   # open julia file in read only mode
    write(ioc, "herr_t save_state (state *s,
   hid_t  file_id)
{
     hid_t group_id;
     herr_t status;

     /* Make Group from simulation time `t` */

     char groupname[50];
     char step_str[10];
     sprintf(step_str, \"%d\", s->step);
     sprintf(groupname, \"%0*d%s\", LEN-(int)strlen(step_str), 0, step_str);

     group_id = H5Gcreate (file_id,
                           groupname,
                           H5P_DEFAULT,
                           H5P_DEFAULT,
                           H5P_DEFAULT);")
    while(!eof(ioj))
        lines = readuntil(ioj, "end\n")
        if startswith(lstrip(lines), "mutable struct")
            empty, lines = split(lines, "mutable struct ")  # strip away the struct
            type, lines... = split(chomp(lines),"\n")   # seperate out the type definition
            # match the type definition to appropriate c types
            clines = @match type begin
                "ReArrays"      =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"status = write_array_dataset(\"\g<var>\", group_id, s->\g<var>, s);\n")
                "halfReArrays"  =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"status = write_Rek_array_dataset(\"\g<var>\", group_id, s->\g<var>, s);\n")
                "dims"          =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"int \g<var>_int = (int)s->\g<var>; \n status = write_int_attribute(\"\g<var>\", group_id, &s->\g<var>);\n")
                "IntParams"     =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"status = write_int_attribute(\"\g<var>\", group_id, &s->\g<var>);\n")
                "doubleParams"  =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"status = write_double_attribute(\"\g<var>\", group_id, &s->\g<var>);\n")
                _ => ""
            end
            write(ioc, join(clines)*"\n")
        end
    end
    write(ioc, "status = H5Gclose (group_id);\nreturn status;\n}\n")
    seekstart(ioj)

    # write the load function
    write(ioc, "state* load_state (hid_t       file_id,
   const char *datafile)
{
    state* s;
    hid_t group_id;
    herr_t status;

    group_id = H5Gopen2 (file_id, datafile, H5P_DEFAULT);")

    # first we need to get the dimension attributes so that we can define the state
    while(!eof(ioj))
        lines = readuntil(ioj, "end\n")
        if startswith(lstrip(lines), "mutable struct dims")
            empty, lines = split(lines, "mutable struct ")  # strip away the struct
            type, lines... = split(chomp(lines),"\n")   # seperate out the type definition
            clines = replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"int \g<var>_int;\n read_int_attribute (\"\g<var>\", group_id, &\g<var>_int);\n ptrdiff_t \g<var> = (ptrdiff_t)\g<var>_int;\n")
            write(ioc, join(clines))
            global dimargs = join(lstrip.(lines) .* ", ")
            global dimargs = rsplit(dimargs,", ",limit=2, keepempty=true)
            global dimargs = join(dimargs)
        end
    end

    seekstart(ioj)
    write(ioc, "s = create_state($dimargs);\n")
    write(ioc, "assert (s != NULL);")
    while(!eof(ioj))
        lines = readuntil(ioj, "end\n")
        if startswith(lstrip(lines), "mutable struct")
            empty, lines = split(lines, "mutable struct ")  # strip away the struct
            type, lines... = split(chomp(lines),"\n")   # seperate out the type definition
            # match the type definition to appropriate c types
            clines = @match type begin
                "ReArrays"      =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"status = read_array_dataset(\"\g<var>\", group_id, &s->\g<var>, s);\n")
                "halfReArrays"  =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"status = read_array_dataset(\"\g<var>\", group_id, &s->\g<var>, s);\n")
                "IntParams"     =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"status = read_int_attribute(\"\g<var>\", group_id, &s->\g<var>);\n")
                "doubleParams"  =>  replace.(lstrip.(lines), r"(?<var>(\w|\d)+)" => s"status = read_double_attribute(\"\g<var>\", group_id, &s->\g<var>);\n")
                _ => ""
            end
            write(ioc, join(clines)*"\n")
        end
    end
    write(ioc, "status = H5Gclose (group_id);\n")
    write(ioc, "return s;\n}\n")
    close(ioc)
end


function generate_dynamics(file, dims)
    ioc = open("dynamics.c","w+")    # generate dynamics file
    write(ioc, "// Define dependencies
#include <stdlib.h>
#include <math.h>
#include <fftw3-mpi.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include \"pfc_header.h\"
#define pi 2.0*acos(0.0)
#define third 1.0 / 3.0\n")

    write_functions(file, dims, ioc)
    close(ioc)
end

# Note that for all for loops I'm going to be assuming that
# 1) we are looping over indices which need to be shifted between julia and c
# 2) since we are looping over indices we are looping over integers
function parse_inline_for(str::String, file)
    # check for a comment in this line and split if appropriate
    try
        line, comment = split(str , "#")
    catch
        line = str
    end
    line = strip(line, ['[',']'])
    strvec = split(line, "for ")
    io = open(file, append)
    braces = 0
    for fornum in 2:length(strvec)
        cline = replace(strvec[fornum], r"(?<index>\w+)( in |=| = |.=|=.)(?<first>\w+):(?<last>\w+)" => s"for(int \g<index> = \g<first> - 1; \g<index>++; \g<index> <\g<last>)")
        write(io, lpad(cline*"\n",4*braces))
        braces+=1
        write(io, lpad("{\n",4*braces))
    end
    try
        write(io, lpad("//"*comment,4*braces))
    catch
    end

    # need handling for various types of indices. Going to assume dimensionality
    # is less than or equal to 3 at least...
    before, after = split(strvec[1], "=")   # split at the equal sign to make life easier as far as spaces go.
    try
        split(before, ",")
        before = replace(r"(?<state_var>\w+)\.(?<type>\w+)\.(?<array_var>\w+)[(?<i>\w+)|((?<i>\w+),(?<j>\w+))|((?<i>\w+),(?<j>\w+),(?<k>\w+))]" => s"\g<state_var>->\g<array_var>[]")
    catch
    end
end

function parse_vectorized(line, ioc, dims, braces)
    # split the vectorized line into the left and right hand sides
    left, right = split(line, "=")
    checkpows = false
    if contains(line, ".^")
        # line = replace(line, ".^" => s"^")
        checkpows = true
    end
    if (checkpows)
        line = parse_power(line)
    end
    # get the operator over which we are vectorizing, and define index max based
    # on whether we are in real space or k space
    statevar, type, array, op = rsplit(left, ".")
    @match dims begin
        2 =>    begin
                    indexends = @match type begin
                        "ReArrays"      => ["s->local_n0" "s->Ny"]
                        "CmxArrays"     => ["s->local_n1" "s->Nx"]
                        "halfReArrays"  => ["s->local_n1" "s->Nx"]
                        _               => [1 1]
                    end
                    # open for loops, define index
                    write(ioc, lpad("for(int i = 0; i < "*indexends[1]*"; i++)\n",4*braces)*lpad("{\n",4braces))
                    braces += 1
                    write(ioc, lpad("for(int j = 0; j < "*indexends[2]*"; j++)\n",4*braces)*lpad("{\n",4*braces))
                    braces += 1
                    @match type begin
                        "ReArrays"      => write(ioc, lpad("index = 2 * i * (("*indexends[2]*">>1)+1) + j;\n", 4*braces))
                        "CmxArrays"     => write(ioc, lpad("index = i * "*indexends[2]*" + j;\n", 4*braces))
                        "halfReArrays"  => write(ioc, lpad("index = i * "*indexends[2]*" + j;\n", 4*braces))
                        _               => write(ioc, lpad("index = i * "*indexends[2]*" + j;\n", 4*braces))
                    end                    # body of the calculation goes here
                    @match type begin
                        "ReArrays"      =>  begin
                                                write(ioc, lpad(statevar*"->"*array*"[index] "*op*"=", 4*braces))
                                                right = replace(right, ".*" => s"*", "./" => s"/", ".+" => s"+", ".-" => s"-", ".^"=>s"^", r"\b(?<func>(\w|\d)+)\.\((?<args>(.)+)\)"=>s"\g<func>(\g<args>)")
                                                if (checkpows)
                                                    right = parse_power(right)
                                                end

                                                right = replace(right, r"(?<state_var>(\w|\d)+)\.ReArrays\.(?<array_var>(\w|\d)+)" => s"\g<state_var>->\g<array_var>[index]")
                                                right = replace(right, r"(?<state_var>(\w|\d)+)\.doubleParams\.(?<param_var>(\w|\d)+)" => s"\g<state_var>->\g<param_var>")
                                                right = vector_standard_replace(right)
                                                write(ioc, lpad(right*";\n", 4*braces))
                                            end
                        "CmxArrays"     =>  begin
                                                write(ioc, lpad(statevar*"->"*array*"[index][0] "*op*"=", 4*braces))
                                                right2 = replace(right, ".*" => s"*", "./" => s"/", ".+" => s"+", ".-" => s"-", ".^"=>s"^", r"\b(?<func>(\w|\d)+)\.\((?<args>(.)+)\)"=>s"\g<func>(\g<args>)")
                                                right = right2
                                                if (checkpows)
                                                    right = parse_power(right)
                                                end

                                                right = replace(right, r"(?<state_var>(\w|\d)+)\.CmxArrays\.(?<array_var>(\w|\d)+)" => s"\g<state_var>->\g<array_var>[index][0]")
                                                right = replace(right, r"(?<state_var>(\w|\d)+)\.doubleParams\.(?<param_var>(\w|\d)+)" => s"\g<state_var>->\g<param_var>")
                                                right = vector_standard_replace(right)
                                                write(ioc, lpad(right*";\n", 4*braces))

                                                write(ioc, lpad(statevar*"->"*array*"[index][1] "*op*"=", 4*braces))
                                                # right = replace(right, ".*" => s"*", "./" => s"/", ".+" => s"+", ".-" => s"-", ".^"=>s"^", r"\b(?<func>(\w|\d)+)\.\((?<args>(.)+)\)"=>s"\g<func>(\g<args>)")
                                                right = right2
                                                if (checkpows)
                                                    right = parse_power(right)
                                                end

                                                right = replace(right, r"(?<state_var>(\w|\d)+)\.CmxArrays\.(?<array_var>(\w|\d)+)" => s"\g<state_var>->\g<array_var>[index][1]")
                                                right = replace(right, r"(?<state_var>(\w|\d)+)\.doubleParams\.(?<param_var>(\w|\d)+)" => s"\g<state_var>->\g<param_var>")
                                                right = vector_standard_replace(right)
                                                write(ioc, lpad(right*";\n", 4*braces))

                                            end
                        "halfReArrays"  =>  begin
                                                write(ioc, lpad(statevar*"->"*array*"[index] "*op*"=", 4*braces))
                                                right = replace(right, ".*" => s"*", "./" => s"/", ".+" => s"+", ".-" => s"-", ".^"=>s"^", r"\b(?<func>(\w|\d)+)\.\((?<args>(.)+)\)"=>s"\g<func>(\g<args>)")

                                                if (checkpows)
                                                    right = parse_power(right)
                                                end

                                                right = replace(right, r"(?<state_var>(\w|\d)+)\.halfReArrays\.(?<array_var>(\w|\d)+)" => s"\g<state_var>->\g<array_var>[index]")
                                                right = replace(right, r"(?<state_var>(\w|\d)+)\.doubleParams\.(?<param_var>(\w|\d)+)" => s"\g<state_var>->\g<param_var>")
                                                right = vector_standard_replace(right)
                                                write(ioc, lpad(right*";\n", 4*braces))

                                            end
                        _               =>  begin
                                                write(ioc, lpad(statevar*"->"*array*" "*op*"=", 4*braces))
                                                right = replace(right, ".*" => s"*", "./" => s"/", ".+" => s"+", ".-" => s"-", ".^"=>s"^", r"\b(?<func>(\w|\d)+)\.\((?<args>(.)+)\)"=>s"\g<func>(\g<args>)")

                                                if (checkpows)
                                                    right = parse_power(right)
                                                end

                                                right = replace(right, r"(?<state_var>\w+)\.(.+)\.(?<array_var>(\w|\d)+)" => s"\g<state_var>->\g<array_var>")
                                                right = vector_standard_replace(right)
                                                write(ioc, lpad(right*";\n", 4*braces))

                                            end
                    end

                    # close the braces
                    braces -= 1
                    write(ioc, lpad("}\n", 4*braces))
                    braces -= 1
                    write(ioc, lpad("}\n", 4*braces))
                end
        3 =>    begin
                    indexends = @match type begin
                        "ReArrays"      => ["s->local_n0" "s->Ny" "s->Nz"]
                        "CmxArrays"     => ["s->local_n1" "s->Nx" "((s->Nz >> 1) + 1)"]
                        "halfReArrays"  => ["s->local_n1" "s->Nx" "((s->Nz >> 1) + 1)"]
                        _               => [0 0 0]
                    end
                    # open for loops, define index
                    write(ioc, lpad("for(int i = 0; i < "*indexends[1]*"; i++)\n",4*braces)*lpad("{\n",4braces))
                    braces += 1
                    write(ioc, lpad("for(int j = 0; j < "*indexends[2]*"; j++)\n",4*braces)*lpad("{\n",4*braces))
                    braces += 1
                    write(ioc, lpad("for(int m = 0; m < "*indexends[3]*"; m++)\n",4*braces)*lpad("{\n",4*braces))
                    braces += 1
                    @match type begin
                        "ReArrays"      => write(ioc, lpad("index = m + 2*(("*indexends[3]*" >> 1) + 1) * (j + "*indexends[2]*"*i);\n", 4*braces))
                        "CmxArrays"     => write(ioc, lpad("index = m + "*indexends[3]*" * (j + "*indexends[2]*"*i);\n", 4*braces))
                        "halfReArrays"  => write(ioc, lpad("index = m + "*indexends[3]*" * (j + "*indexends[2]*"*i);\n", 4*braces))
                        _               => write(ioc, lpad("index = m + "*indexends[3]*" * (j + "*indexends[2]*"*i);\n", 4*braces))
                    end                    # body of the calculation goes here
                    @match type begin
                        "ReArrays"      =>  begin
                                                write(ioc, lpad(statevar*"->"*array*"[index] "*op*"=", 4*braces))
                                                right = replace(right, ".*" => s"*", "./" => s"/", ".+" => s"+", ".-" => s"-", ".^"=>s"^", r"\b(?<func>(\w|\d)+)\.\((?<args>(.)+)\)"=>s"\g<func>(\g<args>)")

                                                if (checkpows)
                                                    right = parse_power(right)
                                                end

                                                right = replace(right, r"(?<state_var>\w+)\.ReArrays\.(?<array_var>(\w|\d)+)" => s"\g<state_var>->\g<array_var>[index]")
                                                right = replace(right, r"(?<state_var>\w+)\.doubleParams\.(?<param_var>(\w|\d)+)" => s"\g<state_var>->\g<param_var>")
                                                # println(right)
                                                right = vector_standard_replace(right)
                                                # println(right)
                                                write(ioc, lpad(right*";\n", 4*braces))
                                            end
                        "CmxArrays"     =>  begin
                                                write(ioc, lpad(statevar*"->"*array*"[index][0] "*op*"=", 4*braces))
                                                right2 = replace(right, ".*" => s"*", "./" => s"/", ".+" => s"+", ".-" => s"-", ".^"=>s"^", r"\b(?<func>(\w|\d)+)\.\((?<args>(.)+)\)"=>s"\g<func>(\g<args>)")
                                                right = right2
                                                if (checkpows)
                                                    right = parse_power(right)
                                                end

                                                right = replace(right, r"(?<state_var>\w+)\.CmxArrays\.(?<array_var>(\w|\d)+)" => s"\g<state_var>->\g<array_var>[index][0]")
                                                right = replace(right, r"(?<state_var>\w+)\.doubleParams\.(?<param_var>(\w|\d)+)" => s"\g<state_var>->\g<param_var>")
                                                right = vector_standard_replace(right)
                                                write(ioc, lpad(right*";\n", 4*braces))

                                                write(ioc, lpad(statevar*"->"*array*"[index][1] "*op*"=", 4*braces))
                                                # right = replace(right, ".*" => s"*", "./" => s"/", ".+" => s"+", ".-" => s"-", ".^"=>s"^", r"\b(?<func>(\w|\d)+)\.\((?<args>(.)+)\)"=>s"\g<func>(\g<args>)")
                                                right = right2
                                                if (checkpows)
                                                    right = parse_power(right)
                                                end

                                                right = replace(right, r"(?<state_var>\w+)\.CmxArrays\.(?<array_var>(\w|\d)+)" => s"\g<state_var>->\g<array_var>[index][1]")
                                                right = replace(right, r"(?<state_var>\w+)\.doubleParams\.(?<param_var>(\w|\d)+)" => s"\g<state_var>->\g<param_var>")
                                                # println(right)
                                                right = vector_standard_replace(right)
                                                # println(right)
                                                write(ioc, lpad(right*";\n", 4*braces))

                                            end
                        "halfReArrays"  =>  begin
                                                # println(checkpows)
                                                write(ioc, lpad(statevar*"->"*array*"[index] "*op*"=", 4*braces))
                                                right = replace(right, ".*" => s"*", "./" => s"/", ".+" => s"+", ".-" => s"-", ".^"=>s"^")

                                                if (checkpows)
                                                    println(right)
                                                    right = parse_power(right)
                                                    println(right)
                                                end

                                                right = replace(right, r"(?<state_var>\w+)\.halfReArrays\.(?<array_var>(\w|\d)+)" => s"\g<state_var>->\g<array_var>[index]")
                                                right = replace(right, r"(?<state_var>\w+)\.doubleParams\.(?<param_var>(\w|\d)+)" => s"\g<state_var>->\g<param_var>")
                                                right = vector_standard_replace(right)
                                                write(ioc, lpad(right*";\n", 4*braces))

                                            end
                        _               =>  begin
                                                write(ioc, lpad(statevar*"->"*array*" "*op*"=", 4*braces))
                                                right = replace(right, ".*" => s"*", "./" => s"/", ".+" => s"+", ".-" => s"-", ".^"=>"^", r"\b(?<func>(\w|\d)+)\.\((?<args>(.)+)\)"=>s"\g<func>(\g<args>)")

                                                if (checkpows)
                                                    right = parse_power(right)
                                                end

                                                right = replace(right, r"(?<state_var>\w+)\.(.+)\.(?<array_var>(\w|\d)+)" => s"\g<state_var>->\g<array_var>")
                                                right = vector_standard_replace(right)
                                                write(ioc, lpad(right*";\n", 4*braces))

                                            end
                    end

                    # close the braces
                    braces -= 1
                    write(ioc, lpad("}\n", 4*braces))
                    braces -= 1
                    write(ioc, lpad("}\n", 4*braces))
                    braces -= 1
                    write(ioc, lpad("}\n", 4*braces))
                end
        _ => write(ioc, "// failed to translate vectorized notation")
    end
    return nothing
end

function parse_local_function_variables(body, args, ioc, braces)
    # where this seems like it would potentially be hardest is determining if
    # the variable has been previously initialized.
    bodycopy = similar(body)
    opset = ["=", "*", "/", "+", "-", "^", ".", ">", "<"]
    others = similar(body)
    for i in eachindex(bodycopy)
        # println(body[i])
        if !(startswith(body[i],"\n") | startswith(body[i],"#")) || (body[i] == "")
            for j in eachindex(opset)
            # left, others... = rsplit(left, opset[j], keepempty=false)
                bodycopy[i], others... = rsplit(body[i], opset[j], keepempty=false)
            end
            # println(bodycopy[i])
            bodycopy[i], others... = rsplit(body[i], opset[1], keepempty=false)
        end
    end

    # Get all variables that are not part of the function arguments
    # cut out the variables that are part of the argument and the return statement
    for i in eachindex(args)
        bodycopy = replace.(bodycopy, args[i]*r"\.(.+)\.(.+)" => s"")
        bodycopy = replace.(bodycopy, r"(.*)\("*args[i]*r"\)" => s"")
        bodycopy = replace.(bodycopy, args[i] => s"")
    end
    # cut out for loops, if, else, end statements, and returns
    bodycopy = replace.(bodycopy, r"for(.+)" => s"","if(.+)"=>s"",r"else(.+)"=>s"",r"return(.+)" => s"", r"(.*)end(.*)"=>s"")
    bodycopy = replace.(bodycopy, r"\s*"=>s"")
    # only lines that remain should be local variable declarations
    # remove empty lines
    # body = setdiff(lstrip.(body), "")
    bodycopy = setdiff(rstrip.(bodycopy), ",")
    bodycopy = setdiff(bodycopy, [""])
    # println(bodycopy)


    # I'm basically going to assume that anything that isn't a index variable is a double

    for i in eachindex(bodycopy)
        if (length(bodycopy[i]) > 0)
            # left, right... = split(body, "=")
            # this is largely for safeties sake, but I'm going to split to check for operators

            cline = lpad("double ".*strip(bodycopy[i]).*";\n",4*braces)
            write(ioc, cline)
        end
    end
    return nothing
end

function standard_replace(line)
    # println(line)
    cline = replace(line,r"(?<statevar>(\w|\d)+)\.(?<type>(ReArrays|CmxArrays|halfReArrays))\.(?<array>(\w|\d)+)\[(.+?)\]" => s"\g<statevar>->\g<array>[index]")
    cline = replace(cline,r"(?<statevar>(\w|\d)+)\.(?<type>(doubleParams|IntParams))\.(?<param>(\w|\d))"=>s"\g<statevar>->\g<param>")
    cline = replace(cline,r"(?<statevar>(\w|\d)+)\.dims\.(?<param>(\w|\d))"=>s"\g<statevar>->\g<param>")
    cline = replace(cline,"2π"=>"2*pi")
    cline = replace(cline,r"end\s*\#(.*)\n"=>s"")
    # println(cline)
    return cline
end
function vector_standard_replace(line)
    # replace vectorized functions to remove the . before the bracket
    test = collect(eachmatch(r"(.*?)\.\(", line))
    for i in eachindex(test)
       line = replace(line, test[i].match=>(test[i].captures)[1]*"(")
   end
    cline = replace(line,r"(?<statevar>(\w|\d)+)\.(?<type>(ReArrays|CmxArrays|halfReArrays))\.(?<array>(\w|\d)+)" => s"\g<statevar>->\g<array>[index]")
    cline = replace(cline,r"(?<statevar>(\w|\d)+)\.(?<type>(doubleParams|IntParams))\.(?<param>(\w|\d)+)"=>s"\g<statevar>->\g<param>")
    cline = replace(cline,"2π"=>"2*pi")
    cline = replace(cline,r"end\s*\#(.*)\n"=>s"")
    return cline
end

function parse_power(line)
    # split up all the powers
    powlines = collect(eachmatch(r".*?\^(\d+)",line))
    # iterate through them
    for m in 1:length(powlines)
        powlines2 = powlines[m].match
        # split off the actual power
        base, power = split(powlines2,r"\s*\.*\^")
        # split off the last open bracket
        base2 = rsplit(base, "(")
        cbase = ""
        numopen = length(base2)
        numclosed = -Inf
        if numopen > 1
            others, lastopen = rsplit(base, "(", limit=2)
            numclosed = length(split(lastopen,")")) - 2
        end
        if ((numopen > 1) & (numclosed >= 0))
            # others, lastopen = rsplit(base, "(", limit=2)
            # get number of closed brackets, which we can then pair up with opens
            # numclosed = length(split(lastopen,")")) - 2
            for n in numopen:-1:(numopen - numclosed)
                cbase = "("*base2[n]*cbase
            end
        elseif numopen > 1
            opsplit = split(lastopen, r"(\+|\-|\*|/)+")
            cbase = opsplit[end]
        else
            # no open brackets, so we need to look for operators to determine the base
            opsplit = split(base, r"(\+|\-|\*|/)+")
            cbase = opsplit[end]
        end
        cline = "pow($cbase,$power)"
        # then recombine with anything remaining in the line
        line = replace(line, "$cbase^$power"=>cline, count=1)
        # println(line)
    end
    return line
end

function parse_function_body(body, args, dims, ioc)
    braces = 1
    checkfors = false
    checkifs = false
    checkwhiles = false
    checkpows = false
    stepin = 0
    # recombine the body so it can be looked at as a totality when necessary.
    body2 = join(body.*"\n")
    # Check the body for for loops
    if contains(body2, r"for+.*?end"ism) || contains(body2, r"for+.*?]"ism)
        # assume that all for loops require an index definition. Not the worst thing
        # if this is not the case.
        write(ioc, lpad("int index;\n",4*braces))
        checkfors = true
    end
    # Also need to check for vectorized loops
    for i in 1:length(args)
        # if (contains(body2, Regex(args[i]*"ReArrays.(+.*?)=")) || contains(body2, Regex(args[i]*"halfReArrays.(+.*?)=")) || contains(body2, Regex(args[i]*"CmxArrays.(+.*?)=")))
        if !checkfors && (contains(body2, args[i]*r"\.ReArrays\.(.*?)=") || contains(body2, args[i]*r"\.halfReArrays\.(.*?)=") || contains(body2, args[i]*r"\.CmxArrays\.(.*?)="))
            write(ioc, lpad("int index;\n",4*braces))
        end
    end
    # Check the body for if statements
    if contains(body2, r"if+.*?end"ism)
        checkifs = true
    end
    # Check the body for while statements
    if contains(body2, r"while+.*?end"ism)
        checkwhiles = true
    end
    # Check the body for powers, because these require a bit more care in handling
    if contains(body2, "^")
        checkpows = true
    end


    # now iterate through each line in the body and parse
    # becuase python is dumb and influenced julia this has to be done with a while loop.
    global i = 1
    # for i in 1:length(body)
    while i <= length(body)
        # case structure for loops, conditionals, returns, comments, etc
        line = lstrip(body[i])
        # first we check for for loops
        if startswith(line, "end")
            # do nothing
            i+=1
            continue
        elseif (checkfors & startswith(line, "for"))
            # set the step in flag until we get to an end call
            # stepin += 1
            # get substring of the for loop
            loop = match(r"for(.*?)end"ism, join(body[i:end].*"\n")).match
            # split at the linebreaks, strip whitespace
            # println(loop)
            loop2 = lstrip.(split(loop, "\n"))

            # get length of the loop
            looplength = length(loop2)-1
            # println("i = $i \t loop = $looplength")
            # println(loop2)

            # get the thing that is being indexed, and it's type
            type = match(r"(?<statevar>\w+)\.(?<type>\w+)\.(.+)\[(.+)\]", loop)[2]
            # check for if we have any part of the body before the start of the next loop
            @match dims begin
                2 => begin
                        subloopcheck = match(r"for(\s|\()(.*)\n(.*)\n(\s*)for", loop)
                        # use type to define index maxima
                        indexends = @match type begin
                            "ReArrays"      => ["s->local_n0" "s->Ny"]
                            "CmxArrays"     => ["s->local_n1" "s->Nx"]
                            "halfReArrays"  => ["s->local_n1" "s->Nx"]
                            _               => [0 0]
                        end
                        # needtranspose = @match type begin
                        #     "ReArrays"      => false
                        #     "CmxArrays"     => true
                        #     "halfReArrays"  => true
                        #     _               => false
                        # end
                        # open for loops, define index
                        write(ioc, lpad("for(int i = 0; i < "*indexends[1]*"; i++)\n",4*braces)*lpad("{\n",4braces))
                        braces += 1
                        # iterate i through the loop
                        for j in 1:looplength
                            global iloop = i + j
                            # println("i = $i \t i2 = $iloop\n")
                            # println(lstrip(body[iloop]))
                            # check if any part of the body occurs before the start of the next loop
                            global counter = 0
                            if subloopcheck != nothing
                                while !(startswith(lstrip(body[iloop]),"for") || subloopcheck == nothing)

                                    # println(iloop)
                                    counter += 1
                                    # write the single loop lines
                                    # if there's a power within the body as a whole, we check this line for a power
                                    if (checkpows & contains(body[iloop],"^"))
                                        # println(body[iloop])
                                        body[iloop] = parse_power(body[iloop])
                                    end
                                    # and finally we handle our generic replacements
                                    cline = standard_replace(body[iloop])
                                    write(ioc, cline*";\n")
                                    if startswith(lstrip(body[iloop+1]),"for")
                                        subloopcheck = nothing

                                    end
                                    global iloop += 1
                                end
                            elseif startswith(lstrip(body[iloop]),"for")
                                write(ioc, lpad("for(int j = 0; j < "*indexends[2]*"; j++)\n",4*braces)*lpad("{\n",4*braces))
                                braces += 1
                                @match type begin
                                    "ReArrays"      => write(ioc, lpad("index = 2 * i * (("*indexends[2]*">>1)+1) + j;\n", 4*braces))
                                    "CmxArrays"     => write(ioc, lpad("index = i * "*indexends[2]*" + j;\n", 4*braces))
                                    "halfReArrays"  => write(ioc, lpad("index = i * "*indexends[2]*" + j;\n", 4*braces))
                                    _               => write(ioc, lpad("index = i * "*indexends[2]*" + j;\n", 4*braces))
                                end
                            end

                            cline = ""

                            j += counter

                            # and write the actual body of the loop
                            # if there's a power within the body as a whole, we check this line for a power
                            if (checkpows & contains(body[iloop],"^"))
                                body[iloop] = parse_power(body[iloop])
                            end
                            # and finally we handle our generic replacements
                            if !startswith(lstrip(body[iloop]),"for") && !startswith(lstrip(body[iloop]),"end")
                                cline = standard_replace(body[iloop])
                                # having handled replacements, write
                                write(ioc, lpad(cline*";\n",4*braces))
                            end
                        end
                        # close loop
                        write(ioc, lpad("}\n",4*braces))
                        braces -= 1
                        write(ioc, lpad("}\n",4*braces))
                        braces -= 1
                        # println("i = $i \t loop = $looplength\n")
                        global i = iloop
                    end
                3 => begin
                        # subloop1check = match(r"for(\s|\()(.*)\n(.*)\n(\s*)for", loop)
                        global loop1length = length(split(split(loop,"for",keepempty=false)[1],r"\n\s*",keepempty=false))
                        global loop2length = length(split(split(loop,"for",keepempty=false)[2],r"\n\s*",keepempty=false))
                        nosubloops = (loop1length+loop2length == 2)
                        # println(loop1length)
                        # println(loop2length)
                        # use type to define index maxima
                        indexends = @match type begin
                            "ReArrays"      => ["s->local_n0" "s->Ny" "s->Nz"]
                            "CmxArrays"     => ["s->local_n1" "s->Nx" "((s->Nz >> 1) + 1)"]
                            "halfReArrays"  => ["s->local_n1" "s->Nx" "((s->Nz >> 1) + 1)"]
                            _               => [0 0 0]
                        end
                        # open for loops, define index
                        write(ioc, lpad("for(int i = 0; i < "*indexends[1]*"; i++)\n",4*braces)*lpad("{\n",4braces))
                        braces += 1
                        # write(ioc, lpad("for(int j = 0; j < "*indexends[2]*"; j++)\n",4*braces)*lpad("{\n",4*braces))
                        # braces += 1
                        # write(ioc, lpad("for(int m = 0; m < "*indexends[3]*"; m++)\n",4*braces)*lpad("{\n",4*braces))
                        # braces += 1

                        # iterate i through the loop
                        global j = 1;
                        while j <= looplength
                            global iloop = i + j
                            # println("i = $i \t i2 = $iloop\n")
                            # println(lstrip(body[iloop]))
                            # check if any part of the body occurs before the start of the next loop
                            global counter = 0
                            if (loop1length > 1)
                                while !(startswith(lstrip(body[iloop]),"for") || loop1length <= 1)

                                    # println(iloop)
                                    counter += 1
                                    # write the single loop lines
                                    # if there's a power within the body as a whole, we check this line for a power
                                    if (checkpows & contains(body[iloop],"^"))
                                        # println(body[iloop])
                                        body[iloop] = parse_power(body[iloop])
                                    end
                                    # and finally we handle our generic replacements
                                    cline = standard_replace(body[iloop])
                                    write(ioc, cline*";\n")
                                    if startswith(lstrip(body[iloop+1]),"for")
                                        global loop1length = 0

                                    end
                                    global iloop += 1
                                    global loop1length -= 1
                                end
                                # global iloop += 1
                                # println("i = $i \t i2 = $iloop\n")

                                if startswith(lstrip(body[iloop]),"for")
                                    write(ioc, lpad("for(int j = 0; j < "*indexends[2]*"; j++)\n",4*braces)*lpad("{\n",4*braces))
                                    braces += 1
                                    global iloop += 1
                                end
                            end

                            global j += counter

                            global counter = 0
                            # cline = ""

                            if (loop2length > 1)
                                while !(startswith(lstrip(body[iloop]),"for") || loop2length <= 1)

                                    # println(iloop)
                                    counter += 1
                                    # write the single loop lines
                                    # if there's a power within the body as a whole, we check this line for a power
                                    if (checkpows & contains(body[iloop],"^"))
                                        # println(body[iloop])
                                        body[iloop] = parse_power(body[iloop])
                                    end
                                    # and finally we handle our generic replacements
                                    cline = standard_replace(body[iloop])
                                    write(ioc, cline*";\n")
                                    if startswith(lstrip(body[iloop+1]),"for")
                                        global loop2length = 0

                                    end
                                    global iloop += 1
                                    global loop2length -= 1
                                end
                                # global iloop += 1
                                # println("i = $i \t i2 = $iloop\n")

                                if startswith(lstrip(body[iloop]),"for")
                                    write(ioc, lpad("for(int m = 0; m < "*indexends[3]*"; m++)\n",4*braces)*lpad("{\n",4*braces))
                                    braces += 1
                                    @match type begin
                                        "ReArrays"      => write(ioc, lpad("index = m + 2*(("*indexends[3]*" >> 1) + 1) * (j + "*indexends[2]*"*i);\n", 4*braces))
                                        "CmxArrays"     => write(ioc, lpad("index = m + "*indexends[3]*" * (j + "*indexends[2]*"*i);\n", 4*braces))
                                        "halfReArrays"  => write(ioc, lpad("index = m + "*indexends[3]*" * (j + "*indexends[2]*"*i);\n", 4*braces))
                                        _               => write(ioc, lpad("index = m + "*indexends[3]*" * (j + "*indexends[2]*"*i);\n", 4*braces))
                                    end
                                    # global iloop += 1
                                    # global counter += 1
                                end
                            end
                            # cline = ""

                            global j += counter

                            if (startswith(lstrip(body[iloop]),"for") & nosubloops)
                                write(ioc, lpad("for(int j = 0; j < "*indexends[2]*"; j++)\n",4*braces)*lpad("{\n",4*braces))
                                braces += 1
                                global counter += 1   # needed so that the loop jumps an extra index
                                # println(j)
                                write(ioc, lpad("for(int m = 0; m < "*indexends[3]*"; m++)\n",4*braces)*lpad("{\n",4*braces))
                                braces += 1
                                # println(j)
                                @match type begin
                                    "ReArrays"      => write(ioc, lpad("index = m + 2*(("*indexends[3]*" >> 1) + 1) * (j + "*indexends[2]*"*i);\n", 4*braces))
                                    "CmxArrays"     => write(ioc, lpad("index = m + "*indexends[3]*" * (j + "*indexends[2]*"*i);\n", 4*braces))
                                    "halfReArrays"  => write(ioc, lpad("index = m + "*indexends[3]*" * (j + "*indexends[2]*"*i);\n", 4*braces))
                                    _               => write(ioc, lpad("index = m + "*indexends[3]*" * (j + "*indexends[2]*"*i);\n", 4*braces))
                                end
                            end

                            # cline = ""

                            global j += counter

                            # and write the actual body of the loop
                            # if there's a power within the body as a whole, we check this line for a power
                            if (checkpows & contains(body[iloop],"^"))
                                body[iloop] = parse_power(body[iloop])
                            end
                            # and finally we handle our generic replacements
                            if !startswith(lstrip(body[iloop]),"for") && !startswith(lstrip(body[iloop]),"end")
                            # println("i = $i \t i2 = $iloop\n")

                                cline = standard_replace(body[iloop])
                                # having handled replacements, write
                                # println(body[iloop])
                                write(ioc, lpad(cline*";\n",4*braces))
                            end
                            global j += 1
                        end
                        # close loop
                        write(ioc, lpad("}\n",4*braces))
                        braces -= 1
                        write(ioc, lpad("}\n",4*braces))
                        braces -= 1
                        write(ioc, lpad("}\n",4*braces))
                        braces -= 1
                        # println("i = $i \t loop = $looplength\n")
                        global i = iloop
                    end
                _ => error("Dimensionality is assumed to be 2 or 3. Check your code.")
            end
            # println(body[i])
            # if startswith(body[i], "end")
            #     # do nothing
            #     i+=1
            # end
        # next we check for if statements
        elseif (checkifs & startswith(body[i], "if"))
            loop = match(r"if+.*?end"ism, body2).match
            # split at the linebreaks, strip whitespace
            loop2 = lstrip(split(loop, "\n"))
            # get length of the loop
            looplength = length(loop2)
            # start the loop
            write(ioc, lpad(loop2[1]*"\n",4*braces)*lpad("{\n",5*(braces+1)))
            for j in 1:looplength
                iloop = i+j
                # if there's a power within the body as a whole, we check this line for a power
                if (checkpows & contains(body[iloop],"^"))
                    body[iloop] = parse_power(body[iloop])
                end
                # handle generic replacements

                if contains(body[i] ,r"\.(\+*|\-*|\**|/*)=")
                    parse_vectorized(body[i], ioc, dims, braces)
                else
                    cline = standard_replace(body[iloop])
                    write(ioc, lpad(cline*";\n",4*braces))
                end
            end
            global i += looplength
            brace -= 1
            write(ioc, lpad("}\n", braces))
        # and while loops
        elseif (checkwhiles & startswith(body[i], "while"))
            # get substring of the for loop
            loop = match(r"while+.*?end"ism, body2).match
            # split at the linebreaks, strip whitespace
            loop2 = lstrip(split(loop, "\n"))
            # get length of the loop
            looplength = length(loop2)
            # start the loop
            write(ioc, lpad(loop2[1]*"\n",4*braces)*lpad("{\n",5*(braces+1)))
            for j in 1:looplength
                iloop = i+j
                # if there's a power within the body as a whole, we check this line for a power
                if (checkpows & contains(body[iloop],"^"))
                    body[iloop] = parse_power(body[iloop])
                end
                # handle generic replacements

                if contains(body[i] ,r"\.(\+*|\-*|\**|/*)=")
                    parse_vectorized(body[i], ioc, dims, braces)
                else
                    cline = standard_replace(body[iloop])
                    write(ioc, lpad(cline*";\n",4*braces))
                end
            end
            global i += looplength
            brace -= 1
            write(ioc, lpad("}\n", braces))

        # then we check if we have a non comparative equals sign, i.e. an assignment
    elseif (contains(body[i], "=") & !(contains(body[i],"==") | contains(body[i],"<=") | contains(body[i],">=") | contains(body[i],"!=")) )
            # check if we are looking at a vectorized input
            # println(body[i])
            if contains(body[i], r"\.(\+*|\-*|\**|/*)=")
                parse_vectorized(body[i], ioc, dims, braces)
                # cline = standard_replace(body[i])
                # write(ioc, lpad(cline*";\n",4*braces))
            # check if we are doing a transform
            elseif contains(body[i], r"\.transforms\.")
                dir = strip(match(r"transforms\.(.*)\s*\*", body[i]).captures[1])
                cline = @match dir begin
                    "forward" =>    replace(body[i], r"(?<statevar>\w+)\.(?<type1>\w+)\.(?<out>.+)\s*=\s*\1\.transforms\.forward\s*\*\s*\1\.(?<type2>\w+)\.(?<in>.+)"=>s"fftw_mpi_execute_dft_r2c(\g<statevar>->forward, \g<statevar>->\g<in>, \g<statevar>->\g<out>);")
                    "back" =>       replace(body[i], r"(?<statevar>\w+)\.(?<type1>\w+)\.(?<out>.+)\s*=\s*\1\.transforms\.back\s*\*\s*\1\.(?<type2>\w+)\.(?<in>.+)"=>s"fftw_mpi_execute_dft_c2r(\g<statevar>->back, \g<statevar>->\g<in>, \g<statevar>->\g<out>);")
                    _ => "// failed to translate transform"
                end
                write(ioc, lpad(cline*"\n",4*braces))
            # otherwise just write the line with standard subs
            else
                cline = standard_replace(body[i])
                write(ioc, lpad(cline*";\n",4*braces))
            end
        elseif startswith(lstrip(body[i]), r"\S*\(.*\)")
            # void outputing function call
            write(ioc, lpad(body[i]*";\n",4*braces))
        elseif startswith(lstrip(body[i]), "return")
            cline = standard_replace(body[i])
            cline = replace(cline, "nothing"=>"")
            write(ioc, lpad(cline*";\n}\n",4*braces))
        else
            # something wonk is happening. Assume comment.
            write(ioc, "// WARNING: COULD NOT MATCH PARSE LINE\n//EXACT JULIA LINE ADDED AS COMMENT\n")
            write(ioc, "// "*body[i]*"\n")
        end
        global i += 1
    end
    return nothing
end

function translate(file, dims)
    # single function to just call the various generate c functions.
    generate_header(file)
    generate_io(file)
    generate_state(file, dims)
    generate_dynamics(file, dims)
end