verdi calculation list -a -p1  # shows Failes jobs
verdi work list -a -p1

verdi group list -A # shows all groups

%verdi work report 102358
2019-02-18 17:28:16 [266 | REPORT]: [102358|PwBaseWorkChain|run_calculation]:
launching PwCalculation<102402> iteration #1
2019-02-18 19:37:27 [326 | REPORT]:
[102358|PwBaseWorkChain|_handle_unexpected_failure]: calculation failure was not
handled
2019-02-18 19:37:27 [327 | REPORT]:
[102358|PwBaseWorkChain|_handle_unexpected_failure]: failure of
PwCalculation<102402> could not be handled, restarting once more
2019-02-18 19:37:28 [328 | REPORT]: [102358|PwBaseWorkChain|run_calculation]:
launching PwCalculation<102588> iteration #2
2019-02-18 21:47:08 [364 | REPORT]:
[102358|PwBaseWorkChain|_handle_unexpected_failure]: calculation failure was not
handled
2019-02-18 21:47:08 [365 | REPORT]:
[102358|PwBaseWorkChain|_handle_unexpected_failure]: failure of
PwCalculation<102588> could not be handled for the second consecutive time
2019-02-18 21:47:13 [366 | REPORT]: [102358|PwBaseWorkChain|on_terminated]: cleaned
remote folders of calculations: 102402 102588


$verdi calculation logshow FAILING_CALC
$verdi calculation outputcat FAILING_CALC

verdi shell
calc = load_node(102588); 
calc.uuid # 6a6de10a-70b5-4992-bd18-9117cd46bc28
In [8]:     calc.get_abs_path()
Out[8]:
u'/home/glensk/aiida/.aiida/repository-al6xxx_test/repository/node/6a/6d/e10a-70b5-4992-bd18-9117cd46bc28'

!!!! this worked best !!!
calc.out.retrieved.get_abs_path()
Out[17]:
u'/home/glensk/aiida/.aiida/repository-al6xxx_test/repository/node/2f/4f/1952-0821-4de4-a734-ea46e114ab41'
