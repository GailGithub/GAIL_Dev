import java.util.*;
public class MeanMC_CLT_Parallel {

	public static double[] main(long nG, int nt) {
		int nGen = (int)nG;		//default value
		long ll, hl = 0;
        int nThreads = nt; //default value
        long [] startPos = new long [nThreads];
        long [] finishPos = new long [nThreads];
		//long time1 = System.nanoTime();	//measuring time
		double value = 0;
        
        GenerateArrayOfRN[] randNumObj = new GenerateArrayOfRN[nThreads];	//defining the array of threads
		
		for (int i=0; i<nThreads; i++)		//setting up the limits for each thread and initializing them as objects
		{
			if(i==0)
				ll = 0;
			else
				ll = hl+1;
			
			if(i==nThreads-1)
				hl = nGen-1;
			else
				hl = hl + nGen/nThreads;
			
            startPos[i] = ll;
            finishPos[i] = hl;
			randNumObj[i] = new GenerateArrayOfRN(ll, hl);
			randNumObj[i].start();
		}
		
        double[] randNumAll= new double[nGen];
		try {
			for (int i=0; i<nThreads; i++)
			{
				randNumObj[i].join();
			}
         
            for (int i=0; i<nThreads; i++){
               System.arraycopy(randNumObj[i].getArray(), 0, randNumAll,(int) startPos[i],(int) randNumObj[i].getArray().length);
            }
            
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		//long time2 = System.nanoTime();
		//long timeTaken = time2 - time1;  
		//System.out.println("Total time taken was " + timeTaken/1000000.0 + " ms");
        

        return randNumAll;
	}

}
