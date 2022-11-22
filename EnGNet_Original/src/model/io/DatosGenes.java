package model.io;

import java.io.*;
import java.util.*;
import util.ReadFromPropertiesSingleton;

public class DatosGenes {
   private HashSet<Gen> listaGenes = new HashSet();
   private ArrayList<String> nombres = new ArrayList();

   public HashSet<Gen> getListaGenes() {
      return this.listaGenes;
   }

   public void setListaGenes(HashSet<Gen> listaGenes) {
      this.listaGenes = listaGenes;
   }

   public void addGen(Gen g) {
      this.listaGenes.add(g);
   }

   public int tamExperimentos() {
      Iterator<Gen> it = this.getListaGenes().iterator();
      Gen g = (Gen)it.next();
      return g.getExperimentos().size();
   }

   public void leerFicheroDeEntrada(String path) throws FichException {
      String linea = null;
      Gen g = null;

      try {
         FileReader file = new FileReader(path);
         BufferedReader br = new BufferedReader(file);
         String[] lista = null;
         linea = br.readLine();

         for(int i = 0; linea != null; ++i) {
            if (i != 0) {
               lista = linea.split("\t");
               g = new Gen(lista[0]);
               this.nombres.add(lista[0]);

               for(int j = 1; j < lista.length; ++j) {
                  if (!lista[j].isEmpty() && !lista[j].equals(" ") && !Float.isNaN(new Float(lista[j]))) {
                     g.addExperimento(new Float(lista[j]));
                  }
               }

               this.addGen(g);
               //System.out.println("Loading line: " + (i + 1));
            }

            linea = br.readLine();
         }

      } catch (FileNotFoundException var9) {
         throw new FichException(var9.getMessage());
      } catch (IOException var10) {
         throw new FichException(var10.getMessage());
      } catch (Exception var11) {
         var11.printStackTrace();
         throw new FichException(var11.getMessage());
      }
   }

   private void trataDeEntrada(String path) {
      String linea = null;
      Object var5 = null;

      try {
         FileReader file = new FileReader(path);
         FileWriter filew = new FileWriter(path + "Salida.txt");
         PrintWriter pw = new PrintWriter(filew);
         BufferedReader br = new BufferedReader(file);
         String[] lista = null;
         linea = br.readLine();

         for(int i = 0; linea != null; ++i) {
            if (i == 0) {
               pw.println(linea);
            } else {
               lista = linea.split(",");
               pw.print(lista[0].substring(1, lista[0].indexOf(40)));

               for(int j = 1; j < lista.length; ++j) {
                  pw.print("\t" + lista[j]);
               }

               pw.print("\n");
               //System.out.println("Loading line: " + (i + 1));
            }

            linea = br.readLine();
         }

         br.close();
         pw.close();
      } catch (Exception var11) {
         var11.printStackTrace();
      }

   }

   public ArrayList<String> listaNombres() {
      return this.nombres;
   }
}