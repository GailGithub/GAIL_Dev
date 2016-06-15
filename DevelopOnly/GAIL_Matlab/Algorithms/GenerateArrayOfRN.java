import java.util.*;

public class GenerateArrayOfRN extends Thread{

	long ll, hl;
	Random rand = new Random();
	ArrayList<Float> arrayRandoms;	
	
	public GenerateArrayOfRN(long lowerLimit, long higherLimit) {
		this.ll = lowerLimit;
		this.hl = higherLimit;
        arrayRandoms = new ArrayList<Float>();	
	}

	public void run ()
	{    
			for (int j = 0; j<=(int)(hl-ll); j++)
			{
				arrayRandoms.add(rand.nextFloat());
			}
	}
    
    public ArrayList<Float> getArray(){
        return arrayRandoms;
    }
}
