import java.util.LinkedList;

/**
 * Created by author.
 */
public class LimitOrder {

    private int type; // 1 for Buy and -1 for Sell
    private Price priceOpened;// is the price at which the limit order was placed
    private double level; // it the price level of the limit order
    private float volume; // is the volume of the order
    private double delta; // is the distance from the price at which the limit order was opened
    private int dcORos; // 1 in case if the order is placed at the level of the expected DC IE, 2 for OS event.
    public LinkedList<LimitOrder> compensatedOrders; // here we should store links to all limit orders aim to be
    // compensated by this one.

    /**
     */
    public LimitOrder(int type, Price priceOpened, double level, float volume, int dcORos, double delta){
        this.type = type;
        this.priceOpened = priceOpened.clone();
        this.level = level;
        this.volume = volume;
        this.delta = delta;
        this.dcORos = dcORos;
        this.compensatedOrders = new LinkedList<>();
    }


    public LimitOrder clone(){
        return new LimitOrder(type, priceOpened, level, volume, dcORos, delta);
    }

    /**
     * Can be called in case if we need to change the level of the order
     * @param level is the new level
     */
    public void setLevel(double level){
        this.level = level;
    }

    public int getType() {
        return type;
    }

    public double getLevel() {
        return level;
    }

    public float getVolume() {
        return volume;
    }

    public int getDcORos() {
        return dcORos;
    }

    public double getDelta() {
        return delta;
    }

    public void addCompenstedOrder(LimitOrder compensatedOrder){
        compensatedOrders.add(compensatedOrder);
    }

    public void cleanCompensatedList(){
        compensatedOrders = new LinkedList<>();
    }

    public void setCompensatedOrders(LinkedList<LimitOrder> compensatedOrders){
        this.compensatedOrders = compensatedOrders;
        volume = computeCompensatedVolume();
    }

    public float computeCompensatedVolume(){
        float compensatedVolume = 0;
        for (LimitOrder aCompensatedOrder : compensatedOrders){
            compensatedVolume += aCompensatedOrder.getVolume();
        }
        return compensatedVolume;
    }


    public double getRelativePnL(){
        double relativePnL = 0;
        for (LimitOrder aCompensatedOrder : compensatedOrders){
            double absPriceMove = (aCompensatedOrder.getLevel() - level) * type;
            if (absPriceMove < 0){
                System.out.println("Negative price move when Sell? " + absPriceMove);
            }
            relativePnL += absPriceMove / aCompensatedOrder.getLevel() * aCompensatedOrder.getVolume();
        }
        return relativePnL;
    }



}
