cd cosmopc
verdi shell
uuid34  ="fe48fbcf-837e-4c2a-9ed4-c7b1e307ba2c" # struct 34 
uuid2786="a8ed0448-d650-4489-b45d-2a8e03409904" # struct 2786
struct 34: ene is -83.069712 hartree in Danieles version == -2260.4418 eV
struct 34: ene is -651.051324359515888 hartree with aiida == -17716.007 eV
struct 34: energy by ** = -17716.00769155667 eV --> correct
uuid = uuid34
load_node(uuid)
calc = load_node(uuid)
calc.out.CALL.out.retrieved.get_abs_path()
struct 34: energy aiida.out = -1302.10268241 rydberg == -17716.008 eV --> correct

struct 2786: ene is -83.1033763 hartree in Danieles version == 2261.3579 eV
struct 2786: ene is -651.400273400486981 hartree with aiida == -17725.503 eV
struct 2786: energy by **  = -17725.50307796094 eV --> correct
uuid = uuid2786
load_node(uuid)
calc = load_node(uuid)
calc.out.CALL.out.retrieved.get_abs_path()
struct 2786: energy aiida.out = -1302.80058051 rydberg == -17725.503 eV --> correct

** = calc.out.CALL.out.output_parameters.get_dict()

USE calc.out.CALL.out.retrieved.get_abs_path() to get the path
    
uuid = uuid34
load_node(uuid)
calc = load_node(uuid)
calc.out.CALL.out.retrieved.get_abs_path()

calc.out.CALL.out.output_parameters.get_dict()
calc.out.CALL.out.output_array.get_array('forces')  # to get the forces
calc.get_abs_path()   # for path to file
calc.out.CALL.get_abs_path()


# abspath 34 -> /home/glensk/aiida/.aiida/repository-al6xxx_test/repository/node/fe/48/fbcf-837e-4c2a-9ed4-c7b1e307ba2c  # is however empty




# abspath 2786-> /home/glensk/aiida/.aiida/repository-al6xxx_test/repository/node/a8/ed/0448-d650-4489-b45d-2a8e03409904 # is also empty
energy hartree = 27.709196390746477 == 754.00558 eV

calc.out.CALL.get_abs_path()
/home/glensk/aiida/.aiida/repository-al6xxx_test/repository/node/a8/bc/df7b-49cd-4147-9d78-c0405a2b7e2c
