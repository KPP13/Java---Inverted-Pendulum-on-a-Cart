package Validator;

// simple functional interface -> lambda expression to check if initial values are in proper periods (SE8)
@FunctionalInterface
public interface Validator {
    boolean check(double xi);
}