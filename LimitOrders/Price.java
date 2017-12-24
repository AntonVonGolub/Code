/**
 * Created by author.
 * Just simple class for prices.
 */

public class Price {

    private float bid;
    private float ask;
    private long time;

    public Price(){}

    public Price(float bid, float ask, long time){
        this.bid = bid;
        this.ask = ask;
        this.time = time;
    }

    public Price clone(){
        return new Price(bid, ask, time);
    }

    public float getAsk() {
        return ask;
    }

    public float getBid() {
        return bid;
    }

    public long getTime() {
        return time;
    }

    public float getMid(){
        return (bid + ask) / 2.0f;
    }

    public void setAsk(int ask){
        this.ask = ask;
    }

    public void setBid(int bid){
        this.bid = bid;
    }

    public void setTime(long time){
        this.time = time;
    }

    public float getSpread(){
        return ask - bid;
    }

}
