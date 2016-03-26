package PhysicalObject;

// abstract class for physical objects with mass, friction coefficients and state variables (position, velocity)
public abstract class PhysicalObject {
    // constraints
    private final double MIN_POS;
    private final double MAX_POS;
    private final double MIN_VEL;
    private final double MAX_VEL;

    // mass
    private final double m;

    // friction coefficient
    private final double d;

    // position / velocity (state variables of physical object)
    private double xp;
    private double xv;

    // constructor is protected to ensure access from outside package
    protected PhysicalObject(final double[] constr, final double mass, final double d, double[] initial) {
        MIN_POS = constr[0];
        MAX_POS = constr[1];
        MIN_VEL = constr[2];
        MAX_VEL = constr[3];
        m = mass;
        this.d = d;
        xp = initial[0];
        xv = initial[1];
    }

    // set state of this object
    public void setState(double[] newState) {
        xp = newState[0];
        xv = newState[1];
    }

    // get state of this object
    public double[] getState() {
        return new double[]{xp, xv};
    }

    // get position constraints (true -> MIN limit, false -> MAX limit)
    public double getPosLimits(boolean min) {
        return (min ? MIN_POS : MAX_POS);
    }

    // get velocity constraints (true -> MIN limit, false -> MAX limit)
    public double getVelLimits(boolean min) {
        return (min ? MIN_VEL : MAX_VEL);
    }

    // get mass
    public double getM() {
        return m;
    }

    // get friction coefficient
    public double getD() {
        return d;
    }

    // rhs mathematical formula - returns RHS value for particular state of the system
    public abstract double[] rhs(double[] state);
}
