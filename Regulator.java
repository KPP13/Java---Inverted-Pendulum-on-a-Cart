package Regulator;

// simple interface for regulator object
@FunctionalInterface
public interface Regulator {
    double calculateControl();        // calculates control value
}