package medidaMatematica;

import java.util.*;
import jsc.correlation.KendallCorrelation;
import jsc.correlation.SpearmanCorrelation;
import jsc.datastructures.PairedData;
import model.io.Gen;
import util.Constantes;

public class Correlacion implements IMedidaMatematica {
   private String tipo;

   public Correlacion(String tipo) {
      this.tipo = Constantes.SPEARMAN;
      this.tipo = tipo;
   }

   public float relacionGenGen(Gen g1, Gen g2) {
      float valor = 0.0F;
      if (this.tipo.equals(Constantes.SPEARMAN)) {
         SpearmanCorrelation sp = new SpearmanCorrelation(new PairedData(this.obtieneValoresGenes(g1.getExperimentos()), this.obtieneValoresGenes(g2.getExperimentos())));
         valor = (float)sp.getR();
      } else if (this.tipo.equals(Constantes.KENDALL)) {
         KendallCorrelation kc = new KendallCorrelation(new PairedData(this.obtieneValoresGenes(g1.getExperimentos()), this.obtieneValoresGenes(g2.getExperimentos())));
         valor = (float)kc.getR();
      }

      return valor;
   }

   double[] obtieneValoresGenes(ArrayList<Float> datos) {
      double[] exp = new double[datos.size()];

      for(int i = 0; i < datos.size(); ++i) {
         exp[i] = ((Float)datos.get(i)).doubleValue();
      }

      return exp;
   }
}