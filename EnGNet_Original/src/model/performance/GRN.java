package model.performance;

import java.io.*;
import java.util.*;
import model.Arco;
import model.Grafo;

public class GRN {
   private ArrayList<String> nodos = new ArrayList();
   private ArrayList<Arco> Arcos = new ArrayList();
   private static GRN instance = null;

   public GRN() {
   }

   public GRN(Grafo g) {
      this.nodos = g.getNombres();
      this.Arcos = g.getAristas();
   }

   public static GRN getInstance() {
      return instance == null ? new GRN() : instance;
   }

   public ArrayList<String> getNodos() {
      return this.nodos;
   }

   public void setNodos(ArrayList<String> nodos) {
      this.nodos = nodos;
   }

   public ArrayList<Arco> getArcos() {
      return this.Arcos;
   }

   public void setArcos(ArrayList<Arco> Arcos) {
      this.Arcos = Arcos;
   }

   public void addArcos(HashSet<Arco> Arcos) {
      this.Arcos.addAll(Arcos);
   }

   public void addNodos(HashSet<String> nodos) {
      this.nodos.addAll(nodos);
   }

   public void addArco(Arco a) {
      this.Arcos.add(a);
   }

   public void volcarAFichero(String path) throws IOException {
      System.out.println("Escribiendo fichero grafo a Fichero");
      System.out.println("Ruta del fichero: " + path);
      System.out.println("Nodos del grafo: " + this.nodos.size());
      System.out.println("Arcos del grafo: " + this.Arcos.size());
      /*FileWriter filew = new FileWriter(path);
      PrintWriter pw = new PrintWriter(filew);
      Iterator<Arco> it = this.Arcos.iterator();
      int j = 0;
      pw.println("Origen\tDestino\tPeso1\tPeso2\tPeso3");

      while(it.hasNext()) {
         Arco Arco = (Arco)it.next();
         pw.println(Arco.toStringESM());
         ++j;
      }

      pw.close();*/
   }

   public void volcarAFicheroNormal(String path) throws IOException {
      System.out.println("Escribiendo fichero grafo a Fichero");
      System.out.println("Ruta del fichero: " + path);
      System.out.println("Nodos del grafo: " + this.nodos.size());
      System.out.println("Arcos del grafo: " + this.Arcos.size());
      FileWriter filew = new FileWriter(path);
      PrintWriter pw = new PrintWriter(filew);
      Iterator<Arco> it = this.Arcos.iterator();
      int j = 0;
      pw.println("Origen\tDestino\tPeso1\tPeso2\tPeso3");

      while(it.hasNext()) {
         Arco Arco = (Arco)it.next();
         pw.println(Arco.toString());
         ++j;
      }

      pw.close();
   }

   public void grafo2Tgf(String path) throws IOException {
      System.out.println("Escribiendo fichero tgf");
      System.out.println("Ruta del fichero: " + path);
      System.out.println("Nodos del grafo: " + this.nodos.size());
      System.out.println("Arcos del grafo: " + this.Arcos.size());
      TreeMap<String, Integer> aris = new TreeMap();
      FileWriter filew = new FileWriter(path);
      PrintWriter pw = new PrintWriter(filew);
      Iterator<String> it = this.nodos.iterator();

      for(int j = 0; it.hasNext(); ++j) {
         String gname = (String)it.next();
         pw.println(j + " " + gname);
         aris.put(gname, j);
      }

      pw.println("#");
      Iterator itA = this.Arcos.iterator();

      while(itA.hasNext()) {
         Arco a = (Arco)itA.next();
         Integer ori = (Integer)aris.get(a.getInicial());
         Integer dest = (Integer)aris.get(a.getTerminal());
         pw.println(ori + " " + dest);
      }

      filew.close();
   }

   public void addNodo(String nodo) {
      this.getNodos().add(nodo);
   }

   private boolean isPresent(Arco a) {
      Iterator it = this.getArcos().iterator();

      while(it.hasNext()) {
         Arco a2 = (Arco)it.next();
         if (a2.igualesAsociacion(a)) {
            return true;
         }
      }

      return false;
   }

   public void removeDuplicatesEdges() {
      for(int i = 0; i < this.getArcos().size() - 1; ++i) {
         Arco a1 = (Arco)this.getArcos().get(i);

         for(int j = i + 1; j < this.getArcos().size(); ++j) {
            Arco a2 = (Arco)this.getArcos().get(j);
            if (a1.igualesAsociacion(a2)) {
               this.getArcos().remove(j);
               System.out.println("edge " + a2.toString() + " removed");
            }
         }
      }

   }
}