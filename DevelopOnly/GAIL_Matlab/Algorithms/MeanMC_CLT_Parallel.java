import java.util.*;
public class MeanMC_CLT_Parallel {

	public static ArrayList<Float> main(long nG, int nt) {
		long nGen = nG;		//default value
		long ll, hl = 0;
        int nThreads = nt;			//default value
		long time1 = System.nanoTime();	//measuring time
		double value = 0;
        
        GenerateArrayOfRN[] randNumObj = new GenerateArrayOfRN[nThreads];	//defining the array of threads
		
		for (int i=0; i<nThreads; i++)		//setting up the limits for each thread and initializing them as objects
		{
			if(i==0)
				ll = 0;
			else
				ll = hl+1;
			
			if(i==nThreads-1)
				hl = nGen;
			else
				hl = hl + nGen/nThreads;
			
			randNumObj[i] = new GenerateArrayOfRN(ll, hl);
			randNumObj[i].start();
		}
		
        ArrayList<Float> randNumAll= new ArrayList<Float>();
		try {
			for (int i=0; i<nThreads; i++)
			{
				randNumObj[i].join();
			}
         
            for (int i=0; i<nThreads; i++){
                   randNumAll.addAll(randNumObj[i].getArray()); 
            }
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		long time2 = System.nanoTime();
		long timeTaken = time2 - time1;  
		System.out.println("Total time taken was " + timeTaken/1000000.0 + " ms");
        return randNumAll;
	}

}
