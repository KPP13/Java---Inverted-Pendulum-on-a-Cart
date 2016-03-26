package InvPendulumSystem;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Locale;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static java.lang.Math.*;

import PhysicalObject.*;
import Regulator.*;
import Solver.*;
import Validator.*;


// inverted pendulum on a cart system class
public final class InvPendSystem {

    // fields
    private boolean regulatorON;                                            // true if regulator is ON
    private PhysicalObject cart;                                            // cart object
    private PhysicalObject pend;                                            // pendulum object
    private final Regulator lq;                                             // lq regulator
    private final Regulator swingUp;                                        // swing up regulator
    private double control;                                                 // control value
    private final double[] k;                                               // constant k coefficients
    private final double[] e;                                               // constant e coefficients
    private final double[] f;                                               // constant f coefficients
    private static double simTime = 0.0;                                    // simulation time
    private final static int N = 4;                                         // order of the system
    public final static int iter = 10000;                                   // number of iterations (default value)


    InvPendSystem(double[] initialState, double[] limits) {

        regulatorON = true;
        cart = new Cart(Arrays.copyOfRange(limits, 0, 4), 0.548, 0.9, Arrays.copyOfRange(initialState, 0, 2));
        pend = new Pendulum(Arrays.copyOfRange(limits, 4, limits.length), 0.131, 0.004, 0.354, 0.0043, Arrays.copyOfRange(initialState, 2, initialState.length));
        control = 0.0;

        final double mc = cart.getM();
        final double mp = pend.getM();
        final double dc = cart.getD();
        final double dp = pend.getD();
        final double g = 9.81;

        double a1 = ((Pendulum)pend).a1;
        double J1 = ((Pendulum)pend).J1;

        k = new double[]{ ((mc + mp) * (mp * pow(a1, 2) + J1)), -(pow(mp, 2) * pow(a1, 2)) };
        e = new double[]{ (-dc * (mp * pow(a1, 2) + J1)), (dp * mp * a1), (mp * a1 * (mp * pow(a1, 2) + J1)),
                (-g * pow(mp, 2) * pow(a1, 2)) / 2 , (mp * pow(a1, 2) + J1) };
        f = new double[]{ (dc * mp * a1), -dp * (mp + mc), (-pow(mp, 2) * pow(a1, 2)) / 2, ((mp + mc) * mp * a1 * g), (-mp * a1) };


        // linear-quadratic regulator (anonymous implementation)
        lq = new Regulator() {
            private final double[] coeffs = {-10.0, -11.0241, -56.5207, -9.6819};

            @Override
            public double calculateControl() {
                double sum = 0;

                for (int i = 0; i < coeffs.length; i++)
                    sum -= getFullState()[i] * coeffs[i];

                return sum;
            }
        };


        // swing-up regulator (anonymous implementation)
        swingUp = new Regulator() {
            @Override
            public double calculateControl() {
                double[] state = getFullState();
                double a1 = ((Pendulum)pend).a1;
                double m1 = pend.getM();

                // swing up regulator algorithm with energy condition
                if (abs(state[3]) <= 0.5) {
                    return -0.2*signum(state[2] - 0.01);
                }
                else {
                    double EK = (m1 * pow(state[3] * a1, 2) / 2 + m1*a1*g * (cos(state[2]) - 1)) + 0.2;

                    return ((EK >= 0) ? 0 : 0.2*signum(state[3] * (abs(state[2]) - PI/2)));             // energy condition
                }
            }
        };
    }


    // cart class
    private class Cart extends PhysicalObject {
        // constructor
        Cart(final double[] constr, final double mass, final double d, double[] initial) {
            super(constr, mass, d, initial);
        }

        // rhs formula
        @Override
        public double[] rhs(double[] x) {
            double[] calc = new double[2];

            calc[0] = x[1];
            calc[1] = (e[0]*x[1] + e[1]*x[3]*cos(x[2]) + e[2]*sin(x[2])*pow(x[3],2) + e[3]*sin(2*x[2]) + e[4]*control) / (k[0] + k[1]*pow(cos(x[2]), 2));

            return calc;
        }
    }


    // pendulum class
    private class Pendulum extends PhysicalObject {
        private final double a1;
        private final double J1;

        // constructor
        Pendulum(final double[] constr, final double mass, final double d, double a, double J, double[] initial) {
            super(constr, mass, d, initial);
            this.a1 = a;
            this.J1 = J;
        }

        // rhs for pendulum
        @Override
        public double[] rhs(double[] x) {
            double[] calc = new double[2];

            calc[0] = x[3];
            calc[1] = (f[0]*x[1]*cos(x[2]) + f[1]*x[3] + f[2]*sin(2*x[2])*pow(x[3],2) + f[3]*sin(x[2]) + f[4]*cos(x[2])*control) / (k[0] + k[1]*pow(cos(x[2]), 2));

            return calc;
        }

        public double getA1() {
            return a1;
        }

        public double getJ1() {
            return J1;
        }
    }


    // RHS of the full system (cart + pendulum)
    public double[] RHS(double[] state) {
        double[] cRHS = cart.rhs(state);
        double[] pRHS = pend.rhs(state);

        return new double[]{ cRHS[0], cRHS[1], pRHS[0], pRHS[1] };
    }


    // returns full state of the system (x1, x2, x3, x4)
    public double[] getFullState() {
        double[] cs = cart.getState();
        double[] ps = pend.getState();

        return new double[]{ cs[0], cs[1], ps[0], ps[1] };
    }


    // calculates control for one step
    public void calcControl(boolean regON) {
        regulatorON = regON;
        double[] state = getFullState();

        // check if the cart is int proper limits
        if (cart.getPosLimits(true) <= state[0] && state[0] <= cart.getPosLimits(false)) {

            if (regulatorON) {
                if (abs(state[2]) <= PI/5)
                    control = lq.calculateControl();            // LQ
                else
                    control = swingUp.calculateControl();       // SU
            }
            else {
                control = 0.0;
            }
        }
        else {
            control = -signum(state[0]);                        // safety condition
        }


        // safe control values
        if (abs(control) >= 0.5)
            control = 0.5*signum(control);

        // convert to N
        control *= 10;
    }


    // solves one step of mathematical model of the system
    public void solveOneStep(boolean regON) {

        calcControl(regON);                             // calculates control
        double[] state = Solver.solveRK45(this);        // solves
        simTime += Solver.h;                            // changes simulation time
        updateState(state);                             // updates state of the system
    }


    // updates state
    private void updateState(double[] state) {
        cart.setState(Arrays.copyOfRange(state,0,2));
        pend.setState(Arrays.copyOfRange(state,2,state.length));

    }


    // returns String with full
    public String getStateInfo() {
        StringBuilder sb = new StringBuilder(100);
        double[] state = getFullState();

        if (simTime == 0)
            sb.append(String.format("%6s\t%15s\t%15s\t%15s\t%15s\t%12s%n", "t[s]", "x1[m]", "x2[m/s]", "x3[rad]", "x4[rad/s]", "F[N]"));

        sb.append(String.format("%6.2f", simTime).replace(",", "."));

        for (int i=0; i<N; i++)
            sb.append(String.format("\t%+10.8e", state[i]).replace(",", "."));

        sb.append(String.format("\t%+8.6e%n", control).replace(",", "."));

        return sb.toString();
    }




    // main method
    public static void main(String[] args) {

        final double[] limits = {-1.5, 1.5, -2, 2, -3.14, 3.14, -5, 5};

        String[] msg = {"Initial cart position (from %4.2f to %4.2f): ", "Initial cart velocity (from %4.2f to %4.2f): ",
                "Initial pendulum position (from %4.2f to %4.2f): ", "Initial pendulum velocity (from %4.2f to %4.2f): "};

        // add limits to text information
        for (int i=0; i<msg.length; i++)
            msg[i] = String.format(msg[i], limits[2*i], limits[2*i+1]).replace(",", ".");


        Double[] initial = new Double[N];   // null references

        // lambda validator array for each state variable
        Validator[] valArr = new Validator[N];

        // create lambda expressions
        for (int i=0; i<4; i++) {
            // j - effectively final for lambda expression
            int j = i;
            valArr[i] = (xi) -> (limits[2*j] <= xi && xi <= limits[2*j + 1]);
        }

        // . as comma in double
        Scanner sc = new Scanner(System.in).useLocale(Locale.US);


        // scan for initial state values and validate
        int i = 0;
        while (i < N) {
            do {
                System.out.println(msg[i]);

                while (!sc.hasNextDouble()) {
                    sc.next();
                    System.out.println(msg[i]);
                }

                Double tmp = sc.nextDouble();
                if (valArr[i].check(tmp)) {
                    initial[i] = tmp;
                }
            } while (initial[i] == null);

            i++;
        }
        sc.close();

        // initial state array
        double[] stateIni = new double[initial.length];
        for (int j=0; j<initial.length; j++)
            stateIni[j] = (double)initial[j];


        // create InvPendSystem object
        InvPendSystem IPSobject = new InvPendSystem(stateIni, limits);


        // prepare a filename to save data
        File file = new File("output1.dat");
        String extension = ".dat";

        // change the index of the file (outputXX.dat), XX - index
        while (file.exists()) {
            // find number of file with regex
            String name = file.getName();
            String regex = "[0-9]+";
            Pattern pat = Pattern.compile(regex);
            Matcher mat = pat.matcher(name);

            if (mat.find()) {
                int index = Integer.parseInt( mat.group() );
                file = new File( "output" + ++index + extension );
            }
        }


        // file not corrupted flag
        boolean fileOK = false;

        // buffered writer
        BufferedWriter bw;

        try {
            bw = new BufferedWriter(new FileWriter(file), 100);
            fileOK = true;
            System.out.println("Data will be stored in " + file.getName() + "file.%n");
        }
        catch (IOException e) {
            System.out.println("Exception occurred:");
            e.printStackTrace();
            System.out.println("Output file will not be created...");
            bw = null;
        }


        // calculate iterations, show data, save to the file
        for (int j=0; j<InvPendSystem.iter; j++) {

            String info = IPSobject.getStateInfo();     // saves the previous state
            System.out.print(info);                   // prints info
            IPSobject.solveOneStep(true);               // solves one step

            // save data to the file (if exists)
            if (bw != null) {
                try {
                    bw.write(info);
                    bw.flush();
                }
                catch (IOException e) {
                    fileOK = false;
                }
            }
        }


        // close file (if exists)
        try {
            if (bw != null)
                bw.close();
        }
        catch (IOException e) {
            e.printStackTrace();
            fileOK = false;
        }

        // final info
        System.out.println(fileOK ? "DONE!" : "DONE, warning: output file might not have been written correctly...");
    }
}
