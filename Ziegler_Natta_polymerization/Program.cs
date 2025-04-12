using System.Transactions;

namespace Ziegler_Natta_polymerization
{
    internal class Program
    {
        static void Main(string[] args)
        {
            double M = 1;
            double totalIngibitor = 0.1;

            double catalystAdditionRate = 0;
            //catalyst state
            double activeCentersAmount = 0;
            double inactiveCentersAmount = 1;
            double poisonedCentersAmount = 0;
            bool poison = false;

            //polymer chains;

            //rate constants
            /*double k_init = 0.1;
            double k_growth = 0.1;
            double k_inhyb = 0.01;*/

            double k_init = 0.5;
            double k_growth = 0.1;
            double k_inhyb = 0.001;

            int iter = 3000;
            int max_chain_length = 1000;
            List<double> inactiveChainsAmount = new List<double>();
            List<double> activeChainsAmount = new List<double>();
            for (int i = 0; i < max_chain_length; i++)
            {
                inactiveChainsAmount.Add(0);
                activeChainsAmount.Add(0);
            }

            int interStep = iter / 100;
            for (int i=0;i < iter;i++)
            {
                //new centers formation
                double dn = k_init * inactiveCentersAmount * M;
                activeCentersAmount += dn;
                inactiveCentersAmount -= dn;
                activeChainsAmount[0] += dn;
                if (inactiveCentersAmount < 0)
                    inactiveCentersAmount = 0;

                //chain continuation
                int max_length = activeChainsAmount.Count;
                for (int j = 0; j < max_chain_length; j++)
                {
                    //chain growth from j to j+1
                    double dn_j_wasting = k_growth * activeChainsAmount[j] * M;
                    //chain growth from j-1 to j
                    double dn_j_producing = 0;
                    if(j!=0)
                        dn_j_producing = k_growth * activeChainsAmount[j - 1] * M;

                    activeChainsAmount[j] -= dn_j_wasting;
                    /*if (j == max_length - 1)
                    {
                        activeChainsAmount.Add(0);
                        inactiveChainsAmount.Add(0);
                    }*/
                    if (j != max_length - 1)
                        activeChainsAmount[j+1] += dn_j_wasting;

                    activeChainsAmount[j] += dn_j_producing;
                    if (j != 0)
                    {
                        activeChainsAmount[j - 1] -= dn_j_producing;
                        if (activeChainsAmount[j - 1] < 0)
                            activeChainsAmount[j - 1] = 0;
                    }

                    if (activeChainsAmount[j] < 0)
                        activeChainsAmount[j] = 0;
                    
                }

                //inhibition of centers
                for (int j = 0; j < max_chain_length; j++)
                {
                    double dn_deactivated = k_inhyb * totalIngibitor * activeChainsAmount[j];
                    activeChainsAmount[j] -= dn_deactivated;
                    if (activeChainsAmount[j] < 0)
                        activeChainsAmount[j] = 0;
                    
                    inactiveChainsAmount[j] += dn_deactivated;
                    activeCentersAmount -= dn_deactivated;
                    if(activeCentersAmount < 0)
                        activeCentersAmount = 0;
                    if (poison)
                        poisonedCentersAmount += dn_deactivated;
                    else
                        inactiveCentersAmount += dn_deactivated;
                }

                //output
                int mod = i % interStep;
                if (mod != 0)
                    continue;
                Console.Write("Iter: " + i + "; ");
                Console.Write(activeCentersAmount + ";");
                Console.Write(inactiveCentersAmount + ";");
                Console.Write("ACTIVE CHAINS;");
                int step = iter / 100;
                if (step == 0)
                    step = 1;
                for (int j = 0; j < max_chain_length; j+=step)
                {
                    Console.Write(activeChainsAmount[j] + ";");
                }
                Console.Write("TOTAL CHAINS;");
                for (int j = 0; j < max_chain_length; j+=step)
                {
                    Console.Write((activeChainsAmount[j]+inactiveChainsAmount[j]) + ";");
                }
                Console.WriteLine();
            }


            //Final output
            Console.WriteLine();
            Console.WriteLine();
            for (int j = 0; j < max_chain_length; j +=5)
            {
                Console.Write(activeChainsAmount[j] + ";");
            }
            Console.WriteLine();
            for (int j = 0; j < max_chain_length; j +=5)
            {
                Console.Write((activeChainsAmount[j] + inactiveChainsAmount[j]) + ";");
            }
            Console.WriteLine();
            Console.WriteLine("k_init = " + k_init);
            Console.WriteLine("k_growth = " + k_growth);
            Console.WriteLine("k_inhyb = " + k_inhyb);
            Console.WriteLine("Iter = " + iter);
            Console.WriteLine("Poisoned = " + poison);
        }
    }
}