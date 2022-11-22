package core;

import java.io.*;
import java.util.*;
import logica.AlgoritmoKruskal;
import medidaMatematica.Correlacion;
import medidaMatematica.IMedidaMatematica;
import medidaMatematica.NMI;
import model.Arco;
import model.Grafo;
import model.Nodo;
import model.io.DatosGenes;
import model.io.FichException;
import model.io.Gen;
import model.performance.GRN;
import util.Constantes;
import util.ReadFromPropertiesSingleton;

public class EnsembleGRN {
   public static void main(String[] args) {
      long inicio = System.currentTimeMillis();
      DatosGenes dg = new DatosGenes();
      ReadFromPropertiesSingleton instance1 = ReadFromPropertiesSingleton.getInstance();
      Properties prop = instance1.getProp();

      try {
         dg.leerFicheroDeEntrada(Constantes.PATH_ENTRADA + System.getProperty("file.separator") + prop.getProperty("inputFile"));
         System.out.println("File loaded...");
         System.out.println("Calculating Ensemble Network...");
         GRN g = generarGrafoCompleto(dg);
         System.out.println("Network computed");
         g.volcarAFichero(Constantes.PATH_ENTRADA + System.getProperty("file.separator") + "grafoCompleto" + prop.getProperty("KendallThreshold") + prop.getProperty("inputFile"));
         Grafo gc = new Grafo(g);
         System.out.println("Pruning of the network starts...");
         Grafo grafoPredominante = AlgoritmoKruskal.aplicarKruskal(gc);
         System.out.println("Adding final relationships...");
         Grafo redFinal = anadirRelacionesGRN(grafoPredominante, dg, prop);
         redFinal.volcarAFichero(Constantes.PATH_ENTRADA + System.getProperty("file.separator") + "finalNetwork" + prop.getProperty("inputFile") + "(" + prop.getProperty("SpearmanThreshold") + ")" + "(" + prop.getProperty("addingThreshold") + ").txt");
         long fin = System.currentTimeMillis();
         System.out.println("Time : " + (fin - inicio) / 1000L + " seconds");
      } catch (FichException var12) {
         System.out.println("An error was occurred: " + var12.getMessage());
      } catch (IOException var13) {
         System.out.println("An error was occurred: " + var13.getMessage());
      }

   }

   private static GRN generarGrafoCompleto(DatosGenes dg) {
      GRN g = GRN.getInstance();
      ReadFromPropertiesSingleton instance1 = ReadFromPropertiesSingleton.getInstance();
      Properties prop = instance1.getProp();
      IMedidaMatematica medida3 = null;
      IMedidaMatematica medida1 = new NMI();
      IMedidaMatematica medida2 = new Correlacion(Constantes.KENDALL);
      medida3 = new Correlacion(Constantes.SPEARMAN);
      int numTotal = dg.getListaGenes().size();
      ArrayList<Gen> lista = new ArrayList(dg.getListaGenes());

      for(int i = 0; i < numTotal - 1; ++i) {
         Gen gene1 = (Gen)lista.get(i);
         if (i % 1000 == 0) {
            System.out.println("Gene:" + gene1.getNombre() + " number: " + i + " of: " + numTotal);
         }

         for(int j = i + 1; j < numTotal; ++j) {
            Gen gene2 = (Gen)lista.get(j);
            int cont = 0;
            if (!gene1.getNombre().equals(gene2.getNombre())) {
               float valor1 = Math.abs(medida1.relacionGenGen(gene1, gene2));
               float valor2 = Math.abs(medida2.relacionGenGen(gene1, gene2));
               float valor3 = Math.abs(medida3.relacionGenGen(gene1, gene2));
               if (valor2 > new Float(prop.getProperty("KendallThreshold"))) {
                  ++cont;
               }

               if (valor3 > new Float(prop.getProperty("SpearmanThreshold"))) {
                  ++cont;
               }

               if (valor1 >= new Float(prop.getProperty("NMIThreshold"))) {
                  ++cont;
               }

               if (cont >= 2) {
                  g.addArco(new Arco(gene1.getNombre(), gene2.getNombre(), valor1, valor2, valor3));
                  if (!g.getNodos().contains(gene1.getNombre())) {
                     g.addNodo(gene1.getNombre());
                  }

                  if (!g.getNodos().contains(gene2.getNombre())) {
                     g.addNodo(gene2.getNombre());
                  }
               }
            }
         }
      }

      return g;
   }

   private static Grafo anadirRelacionesGRN(Grafo grafoPredominante, DatosGenes dg, Properties prop) {
      GRN g = new GRN(grafoPredominante);
      float umbralCor = new Float(prop.getProperty("addingThreshold"));
      int hubThr = new Integer(prop.getProperty("hubThr"));
      if (hubThr < 0) {
         hubThr = 3;
      }

      IMedidaMatematica medida1 = new NMI();
      IMedidaMatematica medida2 = new Correlacion(Constantes.KENDALL);
      IMedidaMatematica medida3 = new Correlacion(Constantes.SPEARMAN);
      Hashtable<String, Nodo> nodos = grafoPredominante.getNodos();
      int k = 1;
      ArrayList<Gen> lista = new ArrayList(dg.getListaGenes());
      int numTotal = dg.getListaGenes().size();

      for(int i = 0; i < numTotal - 1; ++i) {
         Gen gene1 = (Gen)lista.get(i);
         if (nodos.containsKey(gene1.getNombre()) && ((Nodo)nodos.get(gene1.getNombre())).getEnlacesExistentes() > hubThr) {
            System.out.println("Adding relations to relevant node " + k + " : " + gene1.getNombre());

            for(int j = i + 1; j < numTotal; ++j) {
               Gen gene2 = (Gen)lista.get(j);
               //int cont = false;
               if (!gene1.getNombre().equals(gene2.getNombre())) {
                  float valor1 = Math.abs(medida1.relacionGenGen(gene1, gene2));
                  float valor2 = Math.abs(medida2.relacionGenGen(gene1, gene2));
                  float valor3 = Math.abs(medida3.relacionGenGen(gene1, gene2));
                  float media = (valor1 + valor2 + valor3) / 3.0F;
                  if (media >= umbralCor) {
                     g.addArco(new Arco(gene1.getNombre(), gene2.getNombre(), media));
                     if (!g.getNodos().contains(gene2.getNombre())) {
                        g.addNodo(gene2.getNombre());
                     }
                  }
               }
            }
         }

         ++k;
      }

      return grafoPredominante;
   }
}