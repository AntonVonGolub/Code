/**
 * Created by author.
 */


public class Runner {

    private double deltaUp;
    private double deltaDown;
    private double dStarUp;
    private double dStarDown;
    private int mode; // -1 or +1
    private boolean initialized;
    private double extreme;
    private double reference;
    private double expectedDcLevel, expectedOsLevel;


    public Runner(double deltaUp, double deltaDown, double dStarUp, double dStarDown){
        this.deltaUp = deltaUp;
        this.deltaDown = deltaDown;
        this.dStarUp = dStarUp;
        this.dStarDown = dStarDown;
        this.initialized = false;
        this.mode = 1;
    }


    public int run(Price price){
        if(!initialized){
            initialized = true;
            extreme = reference = price.getMid();
            findExpectedDClevel();
            findExpectedOSlevel();
            return 0;
        }

        if( mode == -1 ){
            if( price.getBid() >= expectedDcLevel){
                mode = 1;
                extreme = reference = price.getBid();
                findExpectedDClevel();
                findExpectedOSlevel();
                return 1;
            }
            if( price.getAsk() < extreme ){
                extreme = price.getAsk();
                findExpectedDClevel();
                if( price.getAsk() < expectedOsLevel ){
                    reference = extreme;
                    findExpectedOSlevel();
                    return -2;
                }
            }
        }else if( mode == 1 ){
            if( price.getAsk() <= expectedDcLevel ){
                mode = -1;
                extreme = reference = price.getAsk();
                findExpectedDClevel();
                findExpectedOSlevel();
                return -1;
            }
            if( price.getBid() > extreme ){
                extreme = price.getBid();
                findExpectedDClevel();
                if( price.getBid() > expectedOsLevel ){
                    reference = extreme;
                    findExpectedOSlevel();
                    return 2;
                }
            }
            
        }
        return 0;
    }


    private void findExpectedDClevel(){
        if (mode == -1){
            expectedDcLevel = Math.exp(Math.log(extreme) + deltaUp);
        } else {
            expectedDcLevel = Math.exp(Math.log(extreme) - deltaDown);
        }
    }
    
    
    private void findExpectedOSlevel(){
        if (mode == -1){
            expectedOsLevel = Math.exp(Math.log(reference) - dStarDown);
        } else {
            expectedOsLevel = Math.exp(Math.log(reference) + dStarUp);
        }
    }


    public double getExpectedDcLevel(){
        return expectedDcLevel;
    }


    public double getExpectedOsLevel(){
        return expectedOsLevel;
    }


    public double getExpectedUpperIE(){
        if (expectedDcLevel > expectedOsLevel){
            return expectedDcLevel;
        } else {
            return expectedOsLevel;
        }
    }


    public double getExpectedLowerIE(){
        if (expectedDcLevel < expectedOsLevel){
            return expectedDcLevel;
        } else {
            return expectedOsLevel;
        }
    }


    public int getMode(){
        return mode;
    }


    public double getDeltaUp() {
        return deltaUp;
    }


    public double getDeltaDown() {
        return deltaDown;
    }


    public double getdStarUp() {
        return dStarUp;
    }


    public double getdStarDown() {
        return dStarDown;
    }


    public int getUpperIEtype(){
        if (expectedDcLevel > expectedOsLevel){
            return 1;
        } else {
            return 2;
        }
    }


    public int getLowerIEtype(){
        if (expectedDcLevel < expectedOsLevel){
            return 1;
        } else {
            return 2;
        }
    }


}

