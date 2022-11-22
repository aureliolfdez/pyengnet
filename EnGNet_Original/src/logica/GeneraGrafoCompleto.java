package logica;

import java.io.*;
import java.util.*;
import medidaMatematica.Correlacion;
import medidaMatematica.IMedidaMatematica;
import medidaMatematica.MutualInformation;
import medidaMatematica.SymmetricalUncertainty;
import model.Grafo;
import model.io.DatosGenes;
import model.io.Gen;
import util.Constantes;

public class GeneraGrafoCompleto {
   private static Grafo g = new Grafo();
   private static GeneraGrafoCompleto instance = null;

   private GeneraGrafoCompleto() {
   }

   public static GeneraGrafoCompleto getInstance() {
      if (instance == null) {
         instance = new GeneraGrafoCompleto();
      }

      return instance;
   }

   public Grafo createGrafo() {
      return g;
   }

   public static Grafo generarGrafoCompleto(int inf, int sup, int hilo, DatosGenes dg, String property) {
      try {
         FileWriter filew = new FileWriter(Constantes.PATH_ENTRADA + System.getProperty("file.separator") + "Log-" + hilo + ".txt");
         PrintWriter pw = new PrintWriter(filew);
         IMedidaMatematica medida = null;
         if (!property.equals(Constantes.KENDALL) && !property.equals(Constantes.SPEARMAN)) {
            if (property.equals(Constantes.MI)) {
               medida = new MutualInformation();
            } else if (property.equals(Constantes.SU)) {
               medida = new SymmetricalUncertainty();
            }
         } else {
            medida = new Correlacion(property);
         }

         int i = 1;
         int numTotal = sup - inf;
         ArrayList list = new ArrayList(dg.getListaGenes());

         while(inf < sup) {
            Gen gene1 = (Gen)list.get(inf);
            System.out.println("Hilo " + hilo + "Tratando nodo:" + gene1.getNombre() + " numero: " + i + " de: " + numTotal);
            pw.println("Hilo " + hilo + " Tratando nodo:" + gene1.getNombre() + " numero: " + i + " de: " + numTotal);

            for(int k = inf + 1; k < list.size(); ++k) {
               Gen gene2 = (Gen)list.get(k);
               if (!gene1.getNombre().equals(gene2.getNombre())) {
                  float valor = Math.abs(((IMedidaMatematica)medida).relacionGenGen(gene1, gene2));
                  synchronized(g) {
                     g.ingresarNodo(gene1.getNombre());
                     g.ingresarNodo(gene2.getNombre());
                     g.adicionarEnlace(gene1.getNombre(), gene2.getNombre(), valor);
                  }
               }
            }

            ++inf;
            ++i;
            pw.flush();
         }

         pw.close();
         return g;
      } catch (IOException var17) {
         var17.printStackTrace();
         return null;
      }
   }
}