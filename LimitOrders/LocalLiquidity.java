import java.io.PrintWriter;

/**
 * Created by author.
 */
public class LocalLiquidity {

    double deltaUp, deltaDown;
    double delta;
    double extreme, deltaStar, reference;
    int type;
    boolean initalized;

    double surp;
    double liq;
    double alpha, alphaWeight;
    double H1, H2;


    public LocalLiquidity(double d, double deltaStar, double alpha){
        type = -1; deltaUp = deltaDown = d; this.deltaStar = deltaStar; delta = d;
        initalized = false;
        this.alpha = alpha;
        alphaWeight = Math.exp(-2.0/(alpha + 1.0));
        computeH1H2exp();

    }


    boolean computeH1H2exp(){
        H1 = -Math.exp(-deltaStar/delta)*Math.log(Math.exp(-deltaStar/delta)) - (1.0 - Math.exp(-deltaStar/delta))*Math.log(1.0 - Math.exp(-deltaStar/delta));
        H2 = Math.exp(-deltaStar/delta)*Math.pow(Math.log(Math.exp(-deltaStar/delta)), 2.0) - (1.0 - Math.exp(-deltaStar/delta))*Math.pow(Math.log(1.0 - Math.exp(-deltaStar/delta)), 2.0) - H1*H1;
        return true;
    }


    /**
     * This method should be called with every new price
     * @param price new price
     * @return boolean
     */
    public double computation(Price price){
        int event = run(price);
        if( event != 0 ){
            surp = alphaWeight * (Math.abs(event) == 1 ? 0.08338161 : 2.525729) + (1.0 - alphaWeight) * surp;
            liq = 1.0 - Tools.CumNorm(Math.sqrt(alpha) * (surp - H1) / Math.sqrt(H2));
        }
        return liq;
    }


    /**
     * This is the local runner of the class. It can be delegated to an external class
     * @param price is just a new price
     * @return 1 and -1 if DC IE, 2 or -2 if OS IE, 0 otherwise.
     */
    private int run(Price price){
        if( price == null )
            return 0;

        if( !initalized ){
            type = -1; initalized = true;
            extreme = reference = price.getMid();
            return 0;
        }

        if( type == -1 ){
            if( Math.log(price.getBid()/extreme) >= deltaUp ){
                type = 1;
                extreme = price.getBid();
                reference = price.getBid();
                return 1;
            }
            if( price.getAsk() < extreme ){
                extreme = price.getAsk();
            }
            if( Math.log(reference/extreme) >= deltaStar ){
                reference = extreme;
                return -2;
            }
        }else if( type == 1 ){
            if( Math.log(price.getAsk()/extreme) <= -deltaDown ){
                type = -1;
                extreme = price.getAsk();
                reference = price.getAsk();
                return -1;
            }
            if( price.getBid() > extreme ){
                extreme = price.getBid();
            }
            if( Math.log(reference/extreme) <= -deltaStar ){
                reference = extreme;
                return 2;
            }
        }
        return 0;
    }

}
