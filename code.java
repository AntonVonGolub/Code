/* -- the cleaned up code (without Routes legacy) will be uploaded ASAP -- */

import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.EnumSet;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Random;
import java.util.Set;

import ch.olsen.products.util.configuration.Configuration;
import ch.olsen.products.util.configuration.DoubleProperty;
import ch.olsen.products.util.configuration.EnumProperty;
import ch.olsen.products.util.configuration.StringProperty;
import ch.olsen.routes.atom.Atom;
import ch.olsen.routes.atom.AtomAbstr;
import ch.olsen.routes.atom.AtomException;
import ch.olsen.routes.atom.AtomInput;
import ch.olsen.routes.atom.AtomInputAbstr;
import ch.olsen.routes.atom.AtomOutput;
import ch.olsen.routes.atom.AtomOutputAbstr;
import ch.olsen.routes.atom.RoutesStep;
import ch.olsen.routes.cell.library.LibraryAutoDeploy;
import ch.olsen.routes.data.ArrayDataElement;
import ch.olsen.routes.data.DataElement;
import ch.olsen.routes.data.DataType;
import ch.olsen.routes.data.DoubleDataElement;
import ch.olsen.routes.data.NullDataElement;
import ch.olsen.routes.framework.RoutesFramework;

@LibraryAutoDeploy(name="LykkInvestmentStrategyCode", desc="Lykke Investment Strategy Code", path="Event")
public class Code extends AtomAbstr implements Atom {
	private static final long serialVersionUID = 1L;
	
	/* -- INPUTS -- */
	public AtomInput feed;
	public AtomInput trigger;
	
	/* -- CONFIGURATION -- */
	CodeAtomConfiguration cfg;
	
	public class Runner{
		public double prevExtreme;
		public long prevExtremeTime;
		
		public double prevDC;
		public long prevDCTime;
		
		public double extreme;
		public long extremeTime;
		
		public double deltaUp;
		public double deltaDown;
		public double deltaStarUp, deltaStarDown;
		public double osL;
		public int type;
		public boolean initalized;
		public double reference;
		
		public String fileName;
		
		public Runner(double threshUp, double threshDown, PriceFeedData price, String file, double dStarUp, double dStarDown){
			prevExtreme = price.elems.mid; prevExtremeTime = price.elems.time;
			prevDC = price.elems.mid; prevDCTime = price.elems.time;
			extreme = price.elems.mid; extremeTime = price.elems.time;
			reference = price.elems.mid;
			
			type = -1; deltaUp = threshUp; deltaDown = threshDown; osL = 0.0; initalized = true;
			fileName = new String(file);
			deltaStarUp = dStarUp; deltaStarDown = dStarDown;
		}
		public Runner(double threshUp, double threshDown, double price, String file, double dStarUp, double dStarDown){
			prevExtreme = price; prevExtremeTime = 0;
			prevDC = price; prevDCTime = 0;
			extreme = price; extremeTime = 0;
			reference = price; 
			deltaStarUp = dStarUp; deltaStarDown = dStarDown;
			
			type = -1; deltaUp = threshUp; deltaDown = threshDown; osL = 0.0; initalized = true;
			fileName = new String(file);
		}
		
		public Runner(double threshUp, double threshDown, String file, double dStarUp, double dStarDown){
			deltaUp = threshUp; deltaDown = threshDown;
			initalized = false;
			fileName = new String(file);
			deltaStarUp = dStarUp; deltaStarDown = dStarDown;
		}
		
		public int run(PriceFeedData price){
			if( price == null )
				return 0;
			
			if( !initalized ){
				type = -1; osL = 0.0; initalized = true;
				prevExtreme = price.elems.mid; prevExtremeTime = price.elems.time;
				prevDC = price.elems.mid; prevDCTime = price.elems.time;
				extreme = price.elems.mid; extremeTime = price.elems.time;
				reference = price.elems.mid;
				
				return 0;
			}
			
			if( type == -1 ){
				if( Math.log(price.elems.bid/extreme) >= deltaUp ){
					prevExtreme = extreme;
					prevExtremeTime = extremeTime;
					type = 1;
					extreme = price.elems.ask; extremeTime = price.elems.time;
					prevDC = price.elems.ask; prevDCTime = price.elems.time;
					reference = price.elems.ask;		
					return 1;
				}
				if( price.elems.ask < extreme ){
					extreme = price.elems.ask;
					extremeTime = price.elems.time;
					osL = -Math.log(extreme/prevDC)/deltaDown;
					
					if( Math.log(extreme/reference) <= -deltaStarUp ){
						reference = extreme;
						return -2;
					}
					return 0;
				}
			}else if( type == 1 ){
				if( Math.log(price.elems.ask/extreme) <= -deltaDown ){
					prevExtreme = extreme; 
					prevExtremeTime = extremeTime;
					type = -1;
					extreme = price.elems.bid; extremeTime = price.elems.time;
					prevDC = price.elems.bid; prevDCTime = price.elems.time;
					reference = price.elems.bid;
					return -1;
				}
				if( price.elems.bid > extreme ){
					extreme = price.elems.bid; 
					extremeTime = price.elems.time;
					osL = Math.log(extreme/prevDC)/deltaUp;
					
					if( Math.log(extreme/reference) >= deltaStarDown ){
						reference = extreme;
						return 2;
					}
					return 0;
				}
			}
			
			return 0;
		}
		
		public int run(double price){
			if( !initalized ){
				type = -1; osL = 0.0; initalized = true;
				prevExtreme = price; prevExtremeTime = 0;
				prevDC = price; prevDCTime = 0;
				extreme = price; extremeTime = 0;
				reference = price;
				return 0;
			}
			
			if( type == -1 ){
				if( price - extreme >= deltaUp ){
					prevExtreme = extreme;
					prevExtremeTime = extremeTime;
					type = 1;
					extreme = price; extremeTime = 0;
					prevDC = price; prevDCTime = 0;
					reference = price;
					osL = 0.0;
					
					return 1;
				}
				if( price < extreme ){
					extreme = price;
					extremeTime = 0;
					osL = -(extreme - prevDC);
					if( extreme - reference <= -deltaStarUp ){
						reference = extreme;
						return -2;
					}
					return 0;
				}
			}else if( type == 1 ){
				if( price - extreme <= -deltaDown ){
					prevExtreme = extreme; prevExtremeTime = extremeTime;
					type = -1;
					extreme = price; extremeTime = 0;
					prevDC = price; prevDCTime = 0;
					reference = price;
					osL = 0.0;
					
					return 1;
				}
				if( price > extreme ){
					extreme = price; extremeTime = 0;
					osL = (extreme -prevDC);
					if( extreme - reference >= deltaStarDown ){
						reference = extreme;
						return 2;
					}
					return 0;
				}
			}
			return 0;
		}
	}
	
	public class Liquidity{
		public class Runner{
			public double prevDC;
			public double extreme;
			
			public double deltaUp;
			public double deltaDown;
			public int type;
			public boolean initalized;
			
			public String fileName;
			
			public Runner(double threshUp, double threshDown, PriceFeedData price, String file){
				prevDC = price.elems.mid; 
				extreme = price.elems.mid; 
				
				type = -1; deltaUp = threshUp; deltaDown = threshDown;initalized = true;
				fileName = new String(file);
			}
			public Runner(double threshUp, double threshDown, double price, String file){
				prevDC = price; extreme = price; 
				
				type = -1; deltaUp = threshUp; deltaDown = threshDown;  initalized = true;
				fileName = new String(file);
			}
			
			public Runner(double threshUp, double threshDown, String file){
				deltaUp = threshUp; deltaDown = threshDown;
				initalized = false;
				fileName = new String(file);
			}
			
			public int run(PriceFeedData price){
				if( price == null )
					return 0;
				
				if( !initalized ){
					type = -1; initalized = true;
					prevDC = price.elems.mid;
					extreme = price.elems.mid;
					return 0;
				}
				
				if( type == -1 ){
					if( Math.log(price.elems.bid/extreme) >= deltaUp ){
						type = 1;
						extreme = price.elems.ask; 
						prevDC = price.elems.ask;	
						return 1;
					}
					if( price.elems.ask < extreme ){
						extreme = price.elems.ask;
						return 0;
					}
				}else if( type == 1 ){
					if( Math.log(price.elems.ask/extreme) <= -deltaDown ){
						type = -1;
						extreme = price.elems.bid;
						prevDC = price.elems.bid; 
						return -1;
					}
					if( price.elems.bid > extreme ){
						extreme = price.elems.bid; 
						return 0;
					}
				}
				return 0;
			}
			
			public int run(double price){
				if( !initalized ){
					type = -1; initalized = true;
					prevDC = price;
					extreme = price; 
					return 0;
				}
				
				if( type == -1 ){
					if( price - extreme >= deltaUp ){
						type = 1;
						extreme = price; 
						prevDC = price;
						return 1;
					}
					if( price < extreme ){
						extreme = price;
						return 0;
					}
				}else if( type == 1 ){
					if( price - extreme <= -deltaDown ){
						type = -1;
						extreme = price; 
						prevDC = price; ;
						return 1;
					}
					if( price > extreme ){
						extreme = price;
						return 0;
					}
				}
				return 0;
			}
		}
		
		public Runner[] runner;
		double[] prevState;
		double surp = 0.0, dSurp = 0.0, uSurp = 0.0;
		double liquidity, liquidityUp, liquidityDown; 
		double liqEMA;
		double upLiq, downLiq, diffLiq, diffRaw;
		double H1 = 0.0, H2 = 0.0;
		double d1 = 0.0, d2 = 0.0;
		double alpha, alphaWeight;
		List<Double> mySurprise, downSurprise, upSurprise;
		
		public Liquidity(){};
		public Liquidity(PriceFeedData price, double delta1, double delta2, int lgt){
			double prob = Math.exp(-1.0);
			H1 = -(prob*Math.log(prob) + (1.0 - prob)*Math.log(1.0 - prob));
			H2 = prob*Math.pow(Math.log(prob), 2.0) + (1.0 - prob)*Math.pow(Math.log(1.0 - prob), 2.0) - H1*H1;
			runner = new Runner[lgt];
			prevState = new double[lgt];
			d1 = delta1; d2 = delta2;
	
			getH1nH2(); //skip computation and assign!
			
			runner = new Runner[lgt];
			prevState = new double[lgt];
			
			for( int i = 0; i < runner.length; ++i ){
				runner[i] = new Runner(0.025/100.0 + 0.05/100.0*(double)i, 0.025/100.0 + 0.05/100.0*(double)i, price, "JustFake");
				runner[i].type = (i%2 == 0 ? 1 : -1);
				prevState[i] = (runner[i].type == 1 ? 1 : 0);
			}
			surp = H1; dSurp = H1; uSurp = H1;
			liquidity = 0.5; 
			liqEMA = 0.5;
			
			mySurprise = new LinkedList<Double>();
			downSurprise = new LinkedList<Double>();
			upSurprise = new LinkedList<Double>();
			for( int i = 0; i < 100; ++i ){
				mySurprise.add(new Double(H1));
				downSurprise.add(new Double(H1));
				upSurprise.add(new Double(H1));
			}
			
			//computeLiquidity();
			
			downLiq = 0.5; 
			upLiq = 0.5; 
			diffLiq = 0.5; 
			diffRaw = 0.0;
			alpha = 2.0/(100.0 + 1.0); 
			alphaWeight = Math.exp(-alpha); 
		}
		
		public void getH1nH2(){
			double H1temp = 0.0; double H2temp = 0.0;
			double price = 0.0; 
			alpha = 2.0/(100.0 + 1.0);
			alphaWeight = Math.exp(-alpha);
			runner = new Runner[runner.length];
			for( int i = 0; i < runner.length; ++i ){
				runner[i] = new Runner(0.025/100.0 + 0.05/100.0*(double)i, 0.025/100.0 + 0.05/100.0*(double)i, price, "JustFake");
				runner[i].type = (i%2 == 0 ? 1 : -1);
				prevState[i] = (runner[i].type == 1 ? 1 : 0);
			}
			
			double total1 = 0.0, total2 = 0.0;
			Random rand = new Random(1);
			double dt = 1.0/Math.sqrt(1000000.0);
			double sigma = 0.25; // 25%
			for( int i = 0; i < 100000000; ++i ){
				price += sigma*dt*rand.nextGaussian();
				for( int j= 0; j < runner.length; ++j ){
					if( Math.abs(runner[j].run(price)) == 1 ){ // this is OK for simulated prices
						double myProbs = getProbs(j);
						total1 = total1*alphaWeight + (1.0 - alphaWeight)*(-Math.log(myProbs));
						total2 = total2*alphaWeight + (1.0 - alphaWeight)*Math.pow(Math.log(myProbs), 2.0);
						
						//H1temp = (H1temp*total + -Math.log(myProbs))/(total + 1.0);
						//H2temp = (H2temp*total + Math.pow(Math.log(myProbs), 2.0))/(total + 1.0); 
						//total += 1.0;
					}
				}
			}
			H1 = total1;
			H2 = total2 - H1*H1;
			System.out.println("H1:" + H1 + " H2:" + H2);
		}
		
		public boolean Trigger(PriceFeedData price){
			// -- update values -- 
			boolean doComp = false;
			for( int i = 0; i < runner.length; ++i ){
				int value = runner[i].run(price);
				if( Math.abs(value) == 1 ){
					//double alpha = 2.0/(100.0 + 1.0);
					double myProbs = getProbs(i);
					surp = surp*alphaWeight + (1.0 - alphaWeight)*(-Math.log(myProbs));
					mySurprise.remove(0); mySurprise.add(new Double(-Math.log(myProbs)));
					if( runner[i].type == -1 ){
						dSurp = dSurp*alphaWeight + (1.0 - alphaWeight)*(-Math.log(myProbs));
						downSurprise.remove(0); downSurprise.add(new Double(-Math.log(myProbs)));
					}else if( runner[i].type == 1 ){
						uSurp = uSurp*alphaWeight + (1.0 - alphaWeight)*(-Math.log(myProbs));
						upSurprise.remove(0); upSurprise.add(new Double(-Math.log(myProbs)));
					}
					doComp = true;
				}
			}
			if( doComp ){
				liqEMA = (1.0 - CumNorm(Math.sqrt(100.0)*(surp - H1)/Math.sqrt(H2)));
				upLiq = (1.0 - CumNorm(Math.sqrt(100.0)*(uSurp - H1)/Math.sqrt(H2)));
				downLiq =  (1.0 - CumNorm(Math.sqrt(100.0)*(dSurp - H1)/Math.sqrt(H2)));
				diffLiq = CumNorm(Math.sqrt(100.0)*(uSurp - dSurp)/Math.sqrt(H2));
				diffRaw = Math.sqrt(100.0)*(uSurp-dSurp)/Math.sqrt(H2);
				//computeLiquidity();
			}
			return doComp;
		}
		
		public double getProbs(int i){
			int where = -1;
			for( int j = 1; j < prevState.length; ++j ){
				if( prevState[j] != prevState[0] ){
					where = j;
					break;
				}
			}
			if( i > 0 && where != i ){
				//System.out.println("This should not happen! " + where);
			}
			prevState[i] = (prevState[i] == 1 ? 0 : 1);
			
			if( where == 1 ){
				if( i > 0 ){
					return Math.exp(-(runner[1].deltaDown - runner[0].deltaDown)/runner[0].deltaDown);
				}else{
					return (1.0 - Math.exp(-(runner[1].deltaDown - runner[0].deltaDown)/runner[0].deltaDown));
				}
			}else if( where > 1 ){
				double numerator = 0.0;
				for( int k = 1; k <= where; ++k ){
					numerator -= (runner[k].deltaDown - runner[k-1].deltaDown)/runner[k-1].deltaDown;
				}
				numerator = Math.exp(numerator);
				double denominator = 0.0;
				for( int k = 1; k <= where - 1; ++k ){
					double secVal = 0.0;
					for( int j  = k+1; j <= where; ++j ){
						secVal -=  (runner[j].deltaDown - runner[j-1].deltaDown)/runner[j-1].deltaDown;
					}
					denominator += (1.0 - Math.exp(-(runner[k].deltaDown - runner[k-1].deltaDown)/runner[k-1].deltaDown))*Math.exp(secVal);
				}
				if( i > 0 ){
					return numerator/(1.0 - denominator);
				}else{
					return (1.0 - numerator/(1.0 - denominator));
				}
			}else{
				return 1.0;
			}
		}
		
		// another implementation of the CNDF for a standard normal: N(0,1)
		double CumNorm(double x){
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
		
		
		public boolean computeLiquidity(long deltaT){
			double surpT = 0.0;
			double downSurp = 0.0, upSurp = 0.0;
			
			for( int i = 0; i < mySurprise.size(); ++i ){
				surpT += mySurprise.get(i).doubleValue();
				downSurp += downSurprise.get(i).doubleValue();
				upSurp += upSurprise.get(i).doubleValue();
			}
			
			liquidity = 1.0 - CumNorm((surpT - H1*mySurprise.size())/Math.sqrt(H2*mySurprise.size()));
			liquidityDown = 1.0 - CumNorm((downSurp - H1*downSurprise.size())/Math.sqrt(H2*downSurprise.size()));
			liquidityUp = 1.0 - CumNorm((upSurp - H1*upSurprise.size())/Math.sqrt(H2*upSurprise.size()));
			
			return true;
		}
	};
	
	public class HelpClass{
		long time;
		double price;
		double liq;
		public HelpClass(){};
		public HelpClass(long t, double p, double l){
			time = t; price = p; liq = l;
		}
	};
	
	public class Prices{
		double bid;
		double ask;
		Prices(){};
		Prices(Prices p){
			bid = p.bid;
			ask = p.ask;
		}
		Prices(double b, double a){
			bid = b;
			ask = a;
		}
	};
	
	public class NbDcCounter{
		List<Long> eventList;
		double delta;
		long timeWindow;
		Runner runner;
		
		NbDcCounter(){};
		NbDcCounter(double d, long tW){
			eventList = new LinkedList<Long>();
			delta = d;
			runner = new Runner(delta, delta, "events", delta, delta);
			timeWindow = tW;
		}
		boolean run(PriceFeedData price){
			if( Math.abs(runner.run(price)) == 1 ){
				eventList.add(new Long(price.elems.time));
			}
			
			if( eventList.size() == 0 )
				return true;
			
			while( eventList.get(0).longValue() < price.elems.time - timeWindow )
				eventList.remove(0);
			
			return true;
		}
	};
	
	public class LocalLiquidity{
		double deltaUp, deltaDown;
		double delta;
		double extreme, dStar, reference;
		int type;
		boolean initalized;
		
		double surp, upSurp, downSurp;
		double liq, upLiq, downLiq;
		double alpha, alphaWeight;
		double H1, H2;
		
		LocalLiquidity(){};
		LocalLiquidity(double d, double dUp, double dDown, double dS, double a){
			type = -1; deltaUp = dUp; deltaDown = dDown; dStar = dS; delta = d;
			initalized = false;
			alpha = a;
			alphaWeight = Math.exp(-2.0/(a + 1.0));
			computeH1H2exp(dS);
		}		
		LocalLiquidity(double d,double dUp, double dDown, PriceFeedData price, double dS, double a){
			deltaUp = dUp; deltaDown = dDown; delta = d;
			type = -1;
			extreme = reference = price.elems.mid;
			dStar = dS;
			initalized = true;
			alpha = a;
			alphaWeight = Math.exp(-2.0/(a + 1.0));
			computeH1H2exp(dS);
		}
		boolean computeH1H2exp(double dS){
			H1 = -Math.exp(-dStar/delta)*Math.log(Math.exp(-dStar/delta)) - (1.0 - Math.exp(-dStar/delta))*Math.log(1.0 - Math.exp(-dStar/delta));
			H2 = Math.exp(-dStar/delta)*Math.pow(Math.log(Math.exp(-dStar/delta)), 2.0) - (1.0 - Math.exp(-dStar/delta))*Math.pow(Math.log(1.0 - Math.exp(-dStar/delta)), 2.0) - H1*H1;
			return true;
		}
		// another implementation of the CNDF for a standard normal: N(0,1)
		double CumNorm(double x){
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
		
		public int run(PriceFeedData price){
			if( price == null )
				return 0;
			
			if( !initalized ){
				type = -1; initalized = true;
				extreme = reference = price.elems.mid;
				return 0;
			}
			
			if( type == -1 ){
				if( Math.log(price.elems.bid/extreme) >= deltaUp ){
					type = 1;
					extreme = price.elems.ask;
					reference = price.elems.ask;
					return 1;
				}
				if( price.elems.ask < extreme ){
					extreme = price.elems.ask;
				}
				if( Math.log(reference/extreme) >= dStar ){
					reference = extreme;
					return 2;
				}
			}else if( type == 1 ){
				if( Math.log(price.elems.ask/extreme) <= -deltaDown ){
					type = -1;
					extreme = price.elems.bid; 
					reference = price.elems.bid;
					return -1;
				}
				if( price.elems.bid > extreme ){
					extreme = price.elems.bid; 
				}
				if( Math.log(reference/extreme) <= -dStar ){
					reference = extreme;
					return -2;
				}
			}
			return 0;
		}
		public boolean computation(PriceFeedData price){
			if( price == null )
				return false;
			
			int event = run(price);
			if( event != 0 ){
				surp = alphaWeight*(Math.abs(event) == 1 ? 0.08338161 : 2.525729) + (1.0 - alphaWeight)*surp;
				
				if( event > 0 ){ // down moves
					downSurp = alphaWeight*(event == 1 ? 0.08338161 : 2.525729) + (1.0 - alphaWeight)*downSurp;
				}else if( event < 0 ){ // up moves
					upSurp = alphaWeight*(event == -1 ? 0.08338161 : 2.525729) + (1.0 - alphaWeight)*upSurp;
				}
				
				liq = 1.0 - CumNorm(Math.sqrt(alpha)*(surp - H1)/Math.sqrt(H2)); 
				upLiq = 1.0 - CumNorm(Math.sqrt(alpha)*(upSurp - H1)/Math.sqrt(H2)); 
				downLiq = 1.0 - CumNorm(Math.sqrt(alpha)*(downSurp - H1)/Math.sqrt(H2)); 
			}
			return true;
		}
	};
	
	public class CoastlineTrader{
		double tP; /* -- total position -- */
		List<Double> prices;
		List<Double> sizes;
		
		double profitTarget;
		double pnl, tempPnl;
		double deltaUp, deltaDown, deltaLiq, deltaOriginal;
		double shrinkFlong, shrinkFshort;
		
		double pnlPerc;
		
		int longShort;
		
		boolean initalized;
		Runner runner;
		Runner[][] runnerG;
		
		double increaseLong, increaseShort;
		
		double lastPrice;
		
		double cashLimit;
		String fxRate;
		
		LocalLiquidity liquidity;
		
		CoastlineTrader(){};
		CoastlineTrader(double dOriginal, double dUp, double dDown, double profitT, String FxRate, int lS){
			prices = new LinkedList<Double>();
			sizes = new LinkedList<Double>();
			tP = 0.0; /* -- total position -- */
			
			profitTarget = cashLimit = profitT;
			pnl = tempPnl = pnlPerc = 0.0;
			deltaOriginal = dOriginal;
			deltaUp = dUp; deltaDown = dDown;
			longShort = lS; // 1 for only longs, -1 for only shorts
			shrinkFlong = shrinkFshort = 1.0;
			increaseLong = increaseShort = 0.0;
			
			fxRate = new String(FxRate);
		}
		
		double computePnl(PriceFeedData price){
			// compute PnL with current price
			return 0.0;
		}
		
		double computePnlLastPrice(){
			// compute PnL with last available price
			return 0.0;
		}
		double getPercPnl(PriceFeedData price){
			// percentage PnL
			return 0.0;
		}
		
		boolean tryToClose(PriceFeedData price){
			// PnL target hit implementation
			return false;
		}
		
		boolean assignCashTarget(){
			// implement
			return true;
		}
		
		boolean runPriceAsymm(PriceFeedData price, double oppositeInv){
			if( !initalized ){
				runner = new Runner(deltaUp, deltaDown, price, fxRate, deltaUp, deltaDown);
				
				runnerG = new Runner[2][2];
				
				runnerG[0][0] = new Runner(0.75*deltaUp, 1.50*deltaDown, price, fxRate, 0.75*deltaUp, 0.75*deltaUp);
				runnerG[0][1] = new Runner(0.50*deltaUp, 2.00*deltaDown, price, fxRate, 0.50*deltaUp, 0.50*deltaUp);
				
				runnerG[1][0] = new Runner(1.50*deltaUp, 0.75*deltaDown, price, fxRate, 0.75*deltaDown, 0.75*deltaDown);
				runnerG[1][1] = new Runner(2.00*deltaUp, 0.50*deltaDown, price, fxRate, 0.50*deltaDown, 0.50*deltaDown);
				
				liquidity = new LocalLiquidity(deltaOriginal, deltaUp, deltaDown, price, deltaOriginal*2.525729, 50.0);
				initalized = true;
			}
			
			if( !liquidity.computation(price) ){
				System.out.println("Didn't compute liquidity!");
			}
			
			if( tryToClose(price) ){ /* -- try to close position -- */
				System.out.println("Close");
				return true;
			}
			
			int event = 0;
			
			double fraction = 1.0;
			double size = (liquidity.liq < 0.5 ? 0.5 : 1.0);
			size = (liquidity.liq < 0.1 ? 0.1 : size);
			
			if( longShort == 1 ){ // long positions only
				event = runner.run(price);
				
				if( 15.0 <= tP && tP < 30.0 ){
					event = runnerG[0][0].run(price);
					runnerG[0][1].run(price);
					fraction = 0.5;
				}else if( tP >= 30.0 ){
					event = runnerG[0][1].run(price);
					runnerG[0][0].run(price);
					fraction = 0.25;
				}else{
					runnerG[0][0].run(price); runnerG[0][1].run(price);
				}
				
				if( event < 0 ){
					if( tP == 0.0 ){ // open long position
						int sign = -runner.type;
						if( Math.abs(oppositeInv) > 15.0 ){
							size = 1.0;
							if( Math.abs(oppositeInv) > 30.0 ){
								size = 1.0;
							}
						}
						double sizeToAdd = sign*size; 
						tP += sizeToAdd;
						sizes.add(new Double(sizeToAdd));
						
						prices.add(new Double(sign == 1 ? price.elems.ask : price.elems.bid));
						assignCashTarget();
						System.out.println("Open long");
						
					}else if( tP > 0.0 ){ // increase long position (buy)
						int sign = -runner.type;
						double sizeToAdd = sign*size*fraction*shrinkFlong;
						if( sizeToAdd < 0.0 ){
							System.out.println("How did this happen! increase position but neg size: " + sizeToAdd);
							sizeToAdd = -sizeToAdd; 
						}
						increaseLong += 1.0;
						tP += sizeToAdd;						
						sizes.add(new Double(sizeToAdd));
						
						prices.add(new Double(sign == 1 ? price.elems.ask : price.elems.bid));
						System.out.println("Cascade");
					}
				}else if( event > 0 &&  tP > 0.0 ){ // possibility to decrease long position only at intrinsic events
					double pricE = (tP > 0.0 ? price.elems.bid : price.elems.ask);
					
					for( int i = 1; i < prices.size(); ++i ){
						double tempP = (tP > 0.0 ? Math.log(pricE/prices.get(i).doubleValue()) : Math.log(prices.get(i).doubleValue()/pricE));
						if( tempP >= (tP > 0.0 ? deltaUp : deltaDown) ){
							double addPnl = (pricE - prices.get(i).doubleValue())*sizes.get(i).doubleValue();
							if( addPnl < 0.0 ){
								System.out.println("Descascade with a loss: " + addPnl);
							}
							tempPnl += addPnl;
							tP -= sizes.get(i).doubleValue();
							sizes.remove(i); prices.remove(i);
							increaseLong += -1.0;
							System.out.println("Decascade");
						}
					}
				}
			}else if( longShort == -1 ){ // short positions only
				event = runner.run(price);
				if( -30.0 < tP && tP < -15.0 ){
					event = runnerG[1][0].run(price);
					runnerG[1][1].run(price);
					fraction = 0.5;
				}else if( tP <= -30.0 ){
					event = runnerG[1][1].run(price);
					runnerG[1][0].run(price);
					fraction = 0.25;
				}else{
					runnerG[1][0].run(price); runnerG[1][1].run(price);
				}
				
				if( event > 0 ){
					if( tP == 0.0 ){ // open short position
						int sign = -runner.type;
						if( Math.abs(oppositeInv) > 15.0 ){
							size = 1.0;
							if( Math.abs(oppositeInv) > 30.0 ){
								size = 1.0;
							}
						}
						double sizeToAdd = sign*size;
						if( sizeToAdd > 0.0 ){
							System.out.println("How did this happen! increase position but pos size: " + sizeToAdd);
							sizeToAdd = -sizeToAdd;
						}
						tP += sizeToAdd;
						sizes.add(new Double(sizeToAdd));
						
						prices.add(new Double(sign == 1 ? price.elems.bid : price.elems.ask));
						System.out.println("Open short");
						assignCashTarget();
					}else if( tP < 0.0 ){
						int sign = -runner.type;
						double sizeToAdd = sign*size*fraction*shrinkFshort;
						if( sizeToAdd > 0.0 ){
							System.out.println("How did this happen! increase position but pos size: " + sizeToAdd);
							sizeToAdd = -sizeToAdd;
						}
						
						tP += sizeToAdd;
						sizes.add(new Double(sizeToAdd));
						increaseShort += 1.0;
						prices.add(new Double(sign == 1 ? price.elems.bid : price.elems.ask));
						System.out.println("Cascade");
					}
				}else if( event < 0.0 && tP < 0.0 ){
					double pricE = (tP > 0.0 ? price.elems.bid : price.elems.ask);
					
					for( int i = 1; i < prices.size(); ++i ){
						double tempP = (tP > 0.0 ? Math.log(pricE/prices.get(i).doubleValue()) : Math.log(prices.get(i).doubleValue()/pricE));
						if( tempP >= (tP > 0.0 ? deltaUp : deltaDown) ){
							double addPnl = (pricE - prices.get(i).doubleValue())*sizes.get(i).doubleValue();
							if( addPnl < 0.0 ){
								System.out.println("Descascade with a loss: " + addPnl);
							}
							tempPnl += (pricE - prices.get(i).doubleValue())*sizes.get(i).doubleValue();
							tP -= sizes.get(i).doubleValue();
							sizes.remove(i); prices.remove(i);
							increaseShort += -1.0;
							System.out.println("Decascade");
						}
					}
				}
			}else{
				System.out.println("Should never happen! " + longShort);
			}
			return true;
		}
	};
	
	public class FXrateTrading{
		CoastlineTrader[] coastTraderLong, coastTraderShort;
		String FXrate;
		Liquidity liquidity;
		double currentTime, oneDay;
		boolean init;
		
		FXrateTrading(){};
		FXrateTrading(String rate, int nbOfCoastTraders, double[] deltas){
			currentTime = 1136073600000.0;
			oneDay = 24.0*60.0*60.0*1000.0;
			FXrate = new String(rate);
			coastTraderLong = new CoastlineTrader[nbOfCoastTraders];
			coastTraderShort = new CoastlineTrader[nbOfCoastTraders];
			
			for( int i = 0; i < coastTraderLong.length; ++i ){
				coastTraderLong[i] = new CoastlineTrader(deltas[i], deltas[i], deltas[i], deltas[i], rate.toString(), 1);
				coastTraderShort[i] =  new CoastlineTrader(deltas[i], deltas[i], deltas[i], deltas[i], rate.toString(), -1);
			}
			init = false;
		};
		
		boolean runTradingAsymm(PriceFeedData price){
			if( !init ){
				init = true;
				//liquidity = new Liquidity(price, 0.05/100.0, 0.1/100.0, 50);
			}
			
			//liquidity.Trigger(price);
			//System.out.println("liqEMA: " + liquidity.liqEMA);
			
			for( int i = 0; i < coastTraderLong.length; ++i ){
				coastTraderLong[i].runPriceAsymm(price, coastTraderShort[i].tP);
				coastTraderShort[i].runPriceAsymm(price, coastTraderLong[i].tP);
			}
			
			if( price.elems.time >= currentTime + oneDay ){
				while( currentTime <= price.elems.time )
					currentTime += oneDay;
				
				printDataAsymm(currentTime);
			}
			return true;
		}
		
		boolean printDataAsymm(double time){
			String sep = new String(System.getProperty("file.separator"));
			String folder = new String(sep + "home" + sep + "agolub" + sep + "workspace" + sep + "InvestmentStrategy" + sep + FXrate.toString() + "DataAsymmLiq.dat");
			FileWriter fw = null;
			
			try{
				double totalPos = 0.0, totalShort = 0.0, totalLong = 0.0; double totalPnl = 0.0; double totalPnlPerc = 0.0;  
				fw = new FileWriter(folder, true);
				double price = -1.0;
				for( int i = 0; i < coastTraderLong.length; ++i ){
					if( i == 0 ){
						price = coastTraderLong[i].lastPrice;
					}
					totalLong += coastTraderLong[i].tP;
					totalShort += coastTraderShort[i].tP;
					totalPos += (coastTraderLong[i].tP + coastTraderShort[i].tP);
					totalPnl += (coastTraderLong[i].pnl + coastTraderLong[i].tempPnl + coastTraderLong[i].computePnlLastPrice()
							+ coastTraderShort[i].pnl + coastTraderShort[i].tempPnl + coastTraderShort[i].computePnlLastPrice());
					totalPnlPerc += (coastTraderLong[i].pnlPerc + (coastTraderLong[i].tempPnl + coastTraderLong[i].computePnlLastPrice())/coastTraderLong[i].cashLimit*coastTraderLong[i].profitTarget
							+ coastTraderShort[i].pnlPerc + (coastTraderShort[i].tempPnl + coastTraderShort[i].computePnlLastPrice())/coastTraderShort[i].cashLimit*coastTraderShort[i].profitTarget);
				}
				//double tempSurpScale = Math.sqrt(50)*(liquidity.surp - liquidity.H1)/Math.sqrt(liquidity.H2); 
				fw.append((long)time + "," + totalPnl + "," + totalPnlPerc + "," + totalPos + "," + totalLong + "," + totalShort + "," + price + "\n");
				fw.close();
			}catch(IOException e){
				System.out.println("Failed opening DC thresh file! " + e.getMessage());
				return false;
			}
			return true;
		};
	};
	
	CodeStatus st = new CodeStatus();
	public static class CodeStatus implements Serializable {
		private static final long serialVersionUID = 1L;

		public FXrateTrading[] trading = null;
		public boolean initilized = false;
	}
	
	
	public Code(RoutesFramework framework){
		super(framework, "Engine");
		
		cfg = new CodeAtomConfiguration();
		
		/* -- handling the feed -- */
		//final DataType priceType = PriceFeedData.priceFeedType;
		feed = new AtomInputAbstr(this,"feed", "Input the price feed"){
			private static final long serialVersionUID = 1L;
			public DataType getType() {
				return ArrayDataElement.Factory.buildType(PriceFeedData.priceFeedType);
			}
			public RoutesStep[] receive_internal(DataElement data) throws AtomException {
				if( data == null){
					System.out.println("Null element for input!");
					return null;
				}
				if( st == null ){
					st = new CodeStatus();
				}
				ArrayDataElement ar = data.toArrayDE();
				
				if( !st.initilized ){
					if( !initilize(ar.getAll()) ){
						throw new AtomException("Dynamic Heat Map initilizer failed!");
					}
					return null;
				}else{
					if( !handleScaling(ar.getAll()) ){
						/* -- updating Scaling failed, throw an exception -- */
						throw new AtomException("Dynamic Heat Map update failed!");
					}
					return null;
				}
				
			}
		};
		
		trigger = new AtomInputAbstr(this, "trigger", "End of Data trigger to print"){
			private static final long serialVersionUID = 1L;
			public DataType getType() throws AtomException {
				return NullDataElement.nullDataType;
			}
			public RoutesStep[] receive_internal(DataElement data)
					throws AtomException {
					if( !getPrintOut() ){
						throw new AtomException("Printing failed!");
					}
					return null;
			}
		};	
	}
	
	public boolean initilize(DataElement data[]){
		String[] ccyList = {"AUD_CAD", "AUD_JPY", "AUD_NZD", "AUD_USD", "CAD_JPY", "CHF_JPY", "EUR_AUD", "EUR_CAD", "EUR_CHF",
				"EUR_GBP", "EUR_JPY", "EUR_NZD", "EUR_USD", "GBP_AUD", "GBP_CAD", "GBP_CHF", "GBP_JPY", "GBP_USD", "NZD_CAD",
				"NZD_JPY", "NZD_USD", "USD_CAD", "USD_CHF", "USD_JPY"};
		
		st.trading = new FXrateTrading[ccyList.length];
		double[] deltaS = {0.25/100.0, 0.5/100.0, 1.0/100.0, 1.5/100.0};
		for( int i = 0; i < st.trading.length; ++i ){
			st.trading[i] = new FXrateTrading(ccyList[i], deltaS.length, deltaS);
		}
		st.initilized = true;
		
		for( DataElement price : data ){
			PriceFeedData p = PriceFeedData.cast(price.toAggregatedDE());
			for( int i = 0; i < st.trading.length; ++i ){
				if( st.trading[i].FXrate.equals(p.elems.instrument.toString()) ){
					if( !st.trading[i].runTradingAsymm(p) ){
						System.out.println("Failed at updating!");
						return false;
					}
					break;
				}
			}
		}
		return true;
	}
	
	public boolean handleScaling(DataElement data[]){
		for( DataElement price : data ){
			PriceFeedData p = PriceFeedData.cast(price.toAggregatedDE());
			for( int i = 0; i < st.trading.length; ++i ){
				if( st.trading[i].FXrate.equals(p.elems.instrument.toString()) ){
					if( !st.trading[i].runTradingAsymm(p) ){
						System.out.println("Failed at updating!");
						return false;
					}
					break;
				}
			}
		}
		return true;
	}
	
	public boolean getPrintOut(){
		/*
		for( int k = 0; k < st.array.length; ++k ){
			String ccyName = new String(st.array[k].ccyName);
			String sep = new String(System.getProperty("file.separator"));
			String folder = new String(sep + "home" + sep + "agolub" + sep + "workspace" + sep + "MarketTrend" + sep + ccyName + "MarketTrend.dat");
			FileWriter fw = null;
			
			try{
				fw = new FileWriter(folder, true);
			
				for( int i = 0; i < st.array[k].dcCounter[0].length; ++i ){
					for( int j = 0; j < st.array[k].dcCounter[0].length; ++j ){
						fw.append(st.array[k].dcCounter[i][j] + ",");
					}
					fw.append("\n");
				}
				fw.close();
			}catch(IOException e){
				System.out.println("Failed opening DC thresh file! " + e.getMessage());
				return false;
			}
		}
		*/
		return true;
	}
	
	
	/* -- REMARK:  -- */
	public final String describe() {
		return "CodeAtom";
	}
	
	@Override
	public void reset() {
		st.initilized = false;
		st.trading = null;
	}
	
	public CodeAtomConfiguration getParameters() {
		return cfg;
	}
	
	public enum TypeSpacing{ LOG, LINEAR, BYPASS };
	
	/* -- Sets the Configuration - only two variables, slow and fast moving average arguments -- */
	public static class CodeAtomConfiguration extends 
        Configuration<CodeAtomConfiguration> {

		private static final long serialVersionUID = 0L; /* <- what is this for? */
		
		/* -- input: start day, end day, bin size -- */
		public DoubleProperty startThreshold;
		public DoubleProperty numberOf;
		public DoubleProperty endThreshold;
		//public DoubleProperty ccyNumber;
		public EnumProperty<TypeSpacing> type;
		public StringProperty directory;
		
		
		public CodeAtomConfiguration() {
			
			/* -- starting threshold, and number of DC threshold -- */
			startThreshold = new DoubleProperty("StartingThreshold", "Start Threshold", 0.001, 0.0, Double.MAX_VALUE, false);
			numberOf = new DoubleProperty("NumberofThresholds", "Number Of Thresholds", 10.0, 0.0, Double.MAX_VALUE, false);
			endThreshold = new DoubleProperty("EndThreshold", "End Threshold", 0.05, 0.0, Double.MAX_VALUE, false);
			type = new EnumProperty<TypeSpacing>("SpacingType", "Spacing Type", TypeSpacing.LOG, false);
			//ccyNumber = new DoubleProperty("Number of currencies","NbOfCcy", 20.0, 0.0, Double.MAX_VALUE, false);
			
			// the thresholds are discarded here, as they are set manually.
			
			/* -- take care of this later on -- */
			directory = new StringProperty("Writing directory","WritingDirectory", "myDir", null, false);
		}
		
		public void clear() {
		}

		public String getDescription() {
			return "CodeAtomConfig";
		}

		public String getName() {
			return "codeAtomConfig";
		}
	}
	
	@Override
	public Serializable getManualStatusData() {
		return st;
	}
	@Override
	public void setManualStateData(Serializable value) {
		st = (CodeStatus)value;
	}
}
