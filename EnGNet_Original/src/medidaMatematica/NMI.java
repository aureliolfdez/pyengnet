package medidaMatematica;

import JavaMI.Entropy;
import JavaMI.MutualInformation;
import model.io.Gen;

public class NMI implements IMedidaMatematica {
   public float relacionGenGen(Gen g1, Gen g2) {
      int size = g1.getExperimentos().size();
      float value = 0.0F;
      double[] ge1 = new double[size];
      double[] ge2 = new double[size];

      for(int i = 0; i < size; ++i) {
         ge1[i] = (double)(Float)g1.getExperimentos().get(i);
         ge2[i] = (double)(Float)g2.getExperimentos().get(i);
      }

      try {
         value = 2.0F * (float)MutualInformation.calculateMutualInformation(ge1, ge2) / ((float)Entropy.calculateEntropy(ge1) + (float)Entropy.calculateEntropy(ge2));
      } catch (Exception var8) {
         System.out.println("Error con el gen:  " + g1.getNombre() + " y el gen :" + g2.getNombre());
         value = 0.0F;
      }

      return value;
   }
}