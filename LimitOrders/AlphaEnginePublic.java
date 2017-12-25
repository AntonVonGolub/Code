
public class AlphaEnginePublic {

    public static void main(String[] args) {
        AlphaEngine alphaEngine = new AlphaEngine();
        for (Price price : priceFeed){ // please use your own price feed
            alphaEngine.run(price);
        }
    }

}