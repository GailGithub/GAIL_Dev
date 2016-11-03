import java.util.*;

public class GenerateArrayOfRN extends Thread{

	long ll, hl;
	Random rand = new Random();
	double[] arrayRandoms;	
	
	public GenerateArrayOfRN(long lowerLimit, long higherLimit) {
		this.ll = lowerLimit;
		this.hl = higherLimit;
        arrayRandoms = new double [(int)(hl-ll)+1];	
	}

	public void run ()
	{    
			for (int j = 0; j<=(int)(hl-ll); j++)
			{
				arrayRandoms[j] = rand.nextDouble();
               
			}
	}
    
    public double[] getArray(){
        return arrayRandoms;
    }
}
