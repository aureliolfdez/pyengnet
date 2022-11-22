package logica;

import java.util.*;
import model.Arco;
import model.Enlace;
import model.Grafo;
import model.Nodo;

public class AlgoritmoKruskal {
   public static Grafo aplicarKruskal(Grafo grafo) {
      Grafo arbol = new Grafo();
      ArrayList<String> nodos = grafo.getNombres();

      for(int j = 0; j < nodos.size(); ++j) {
         arbol.ingresarNodo((String)nodos.get(j));
      }

      ArrayList<Arco> L = (ArrayList)grafo.getAristas().clone();
      Arco pro = (Arco)L.get(0);
      arbol.adicionarEnlace(pro.getInicial(), pro.getTerminal(), pro.getPeso());
      L.remove(pro);

      for(int i = 0; L.size() != 0; ++i) {
         if (i % 250 == 0) {
            System.out.println("Working on cicle " + i + " of " + L.size());
         }

         pro = (Arco)L.get(0);
         if (!HayCiclo(arbol, pro, arbol.getNodo(pro.getTerminal()), pro.getTerminal())) {
            arbol.adicionarEnlace(pro.getInicial(), pro.getTerminal(), pro.getPeso());
         }

         L.remove(pro);
      }

      return arbol;
   }

   public static boolean HayCiclo(Grafo g, Arco aVerificar, Nodo terminal, String N) {
      ArrayList<Enlace> aux = terminal.getEnlaces();
      if (aux.size() == 0) {
         return false;
      } else if (terminal.existeEnlace(aVerificar.getInicial()) != -1) {
         return true;
      } else {
         for(int i = 0; i < aux.size(); ++i) {
            Enlace nodo = (Enlace)aux.get(i);
            if (!nodo.getDestino().equals(N) && HayCiclo(g, aVerificar, g.getNodo(nodo.getDestino()), terminal.getNombre())) {
               return true;
            }
         }

         return false;
      }
   }
}
