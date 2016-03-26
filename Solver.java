package Solver;
import InvPendulumSystem.InvPendSystem;


// RK45 solver for time invariant system x' = f(x) (modified algorithm)
public final class Solver {
    // step
    public static double h = 0.01;

    // RK45 solver
    public static double[] solveRK45(InvPendSystem obj) {

        // actual state
        double[] state = obj.getFullState();

        // number of equations
        int n = state.length;

        // k coefficients
        double[][] k = new double[4][];

        // temporary state
        double[] new_state = new double[n];

        // RHS computation result in each step of the algorithm
        double[] tmp_res;


        // algorithm starts here
        for (int i=0; i<k.length; i++) {

            k[i] = new double[n];
            tmp_res = obj.RHS( i==0 ? state : new_state );

            for (int j=0; j<n; j++) {

                k[i][j] = h * tmp_res[j];

                if (i != k.length-1) {
                    new_state[j] = state[j] + 0.5 * k[i][j]*h;
                    continue;
                }

                new_state[j] = (k[0][j] + 2*k[1][j] + 2*k[2][j] + k[3][j]) / 6;
                state[j] += new_state[j];
            }
        }

        return state;
    }
}
