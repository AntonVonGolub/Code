
/**
 * Created by author.
 * The class Tools holds a set of auxiliary functions which can be used in order to simplify work on the general ideas.
 */
public class Tools {


    /**
     * The method computes cumulative value of the normal distribution with parameter x
     * @param x is the x coordinate of the normal distribution
     * @return sum of the cumulative normal distribution
     */
    public static double CumNorm(double x){
        // protect against overflow
        if (x > 6.0)
            return 1.0;
        if (x < -6.0)
            return 0.0;
        double b1 = 0.31938153;
        double b2 = -0.356563782;
        double b3 = 1.781477937;
        double b4 = -1.821255978;
        double b5 = 1.330274429;
        double p = 0.2316419;
        double c2 = 0.3989423;
        double a = Math.abs(x);
        double t = 1.0 / (1.0 + a * p);
        double b = c2*Math.exp((-x)*(x/2.0));
        double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
        n = 1.0-b*n;
        if ( x < 0.0 )
            n = 1.0 - n;
        return n;
    }


}
