package model;

import java.io.*;
import java.util.*;
import model.performance.GRN;

public class Grafo {
   ArrayList<String> nombres;
   ArrayList<Arco> aristas;
   Hashtable<String, Nodo> nodos;

   public Grafo() {
      this.nombres = new ArrayList();
      this.nodos = new Hashtable();
      this.aristas = new ArrayList();
   }

   public Grafo(GRN grn) {
      this.nombres = new ArrayList(grn.getNodos());
      this.aristas = grn.getArcos();
      this.nodos = new Hashtable();
      Iterator it = grn.getNodos().iterator();

      while(it.hasNext()) {
         this.ingresarNodo((String)it.next());
      }

   }

   public void ingresarNodo(String nombre) {
      if (this.nodos.get(nombre) == null) {
         this.nombres.add(nombre);
         this.nodos.put(nombre, new Nodo(nombre));
      }

   }

   public void adicionarEnlace(String nodoInicial, String nodoTerminal, float peso) {
      Arco nuevo = new Arco(nodoInicial, nodoTerminal, peso);
      int i = this.buscarIndice(nuevo.getPeso());
      if (i == -1) {
         this.aristas.add(nuevo);
      } else {
         this.aristas.add(i, nuevo);
      }

      ((Nodo)this.nodos.get(nodoInicial)).agregarEnlace(nodoTerminal, (double)peso);
      ((Nodo)this.nodos.get(nodoTerminal)).agregarEnlace(nodoInicial, (double)peso);
   }

   public boolean busarArista(Arco arco) {
      for(int i = 0; i < this.aristas.size(); ++i) {
         Arco otro = (Arco)this.aristas.get(i);
         if (arco.getInicial().equals(otro.getInicial()) && arco.getTerminal().equals(otro.getTerminal()) && arco.getPeso() == otro.getPeso()) {
            this.aristas.remove(otro);
            return true;
         }
      }

      return false;
   }

   public boolean existeArista(Arco arco) {
      for(int i = 0; i < this.aristas.size(); ++i) {
         Arco otro = (Arco)this.aristas.get(i);
         if (arco.igualesAsociacion(otro)) {
            return true;
         }
      }

      return false;
   }

   public int buscarIndice(float peso) {
      for(int i = 0; i < this.aristas.size(); ++i) {
         if (peso > ((Arco)this.aristas.get(i)).getPeso()) {
            return i;
         }
      }

      return -1;
   }

   public Hashtable getNodos() {
      return this.nodos;
   }

   public void setNodos(Hashtable<String, Nodo> muchos) {
      this.nodos = muchos;
   }

   public ArrayList<String> getNombres() {
      return this.nombres;
   }

   public Nodo getNodo(String nombre) {
      return (Nodo)this.nodos.get(nombre);
   }

   public ArrayList<Arco> getAristas() {
      return this.aristas;
   }

   public void setAristas(ArrayList<Arco> aristas) {
      this.aristas = aristas;
   }

   public void setNombres(ArrayList<String> nombres) {
      this.nombres = nombres;
   }

   public void volcarAFichero(String path) throws IOException {
      System.out.println("Writing network into a file");
      System.out.println("File path: " + path);
      System.out.println("Number of genes of the network: " + this.nodos.size());
      System.out.println("Number of edges of the network: " + this.aristas.size());
      /*FileWriter filew = new FileWriter(path);
      PrintWriter pw = new PrintWriter(filew);
      Iterator<Arco> it = this.aristas.iterator();
      int j = 0;
      pw.println("Gene1\tGene2\tWeight");

      while(it.hasNext()) {
         Arco arista = (Arco)it.next();
         pw.println(arista.toString());
         ++j;
      }

      pw.close();*/
   }
}