#!/bin/sh
diff FLa4a_tmb.Rcheck/tests/addition.Rout FLa4a.Rcheck/tests/addition.Rout | cat > addition.diff
diff FLa4a_tmb.Rcheck/tests/l2a.Rout FLa4a.Rcheck/tests/l2a.Rout | cat > l2a.diff
diff FLa4a_tmb.Rcheck/tests/m.Rout FLa4a.Rcheck/tests/m.Rout | cat > m.diff
diff FLa4a_tmb.Rcheck/tests/predict.Rout FLa4a.Rcheck/tests/predict.Rout | cat > predict.diff
diff FLa4a_tmb.Rcheck/tests/residuals.Rout FLa4a.Rcheck/tests/residuals.Rout | cat > residuals.diff
diff FLa4a_tmb.Rcheck/tests/sca.Rout FLa4a.Rcheck/tests/sca.Rout | cat > sca.diff
diff FLa4a_tmb.Rcheck/tests/scas.Rout FLa4a.Rcheck/tests/scas.Rout | cat > scas.diff
diff FLa4a_tmb.Rcheck/tests/simulate.Rout FLa4a.Rcheck/tests/simulate.Rout | cat > simulate.diff
diff FLa4a_tmb.Rcheck/tests/wcsam.Rout FLa4a.Rcheck/tests/wcsam.Rout | cat > wcsam.diff
diff FLa4a_tmb.Rcheck/00check.log FLa4a.Rcheck/00check.log | cat > 00check.diff


