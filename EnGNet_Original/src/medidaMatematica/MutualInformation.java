package medidaMatematica;

import model.io.Gen;

public class MutualInformation implements IMedidaMatematica {
   public float relacionGenGen(Gen g1, Gen g2) {
      int size = g1.getExperimentos().size();
      double[] ge1 = new double[size];
      double[] ge2 = new double[size];

      for(int i = 0; i < size; ++i) {
         ge1[i] = (double)(Float)g1.getExperimentos().get(i);
         ge2[i] = (double)(Float)g2.getExperimentos().get(i);
      }

      float value = (float)JavaMI.MutualInformation.calculateMutualInformation(ge1, ge2);
      return value;
   }
}