package model;

import java.awt.Component;
import java.util.*;
import javax.swing.JOptionPane;

public class Nodo {
   private String nombre;
   private ArrayList<Enlace> enlaces;
   private int enlacesExistentes;

   public ArrayList<Enlace> getEnlaces() {
      return this.enlaces;
   }

   public Nodo(String newName) {
      this.nombre = newName;
      this.enlacesExistentes = -1;
      this.enlaces = new ArrayList();
   }

   public int getEnlacesExistentes() {
      return this.enlacesExistentes;
   }

   public String getNombre() {
      return this.nombre;
   }

   public void agregarEnlace(String enlazar, double peso) {
      if (this.enlacesExistentes == -1) {
         this.enlaces.add(new Enlace(enlazar, peso));
         ++this.enlacesExistentes;
      } else {
         int posicion = this.existeEnlace(enlazar);
         if (posicion == -1) {
            this.enlaces.add(new Enlace(enlazar, peso));
            ++this.enlacesExistentes;
         }
      }

   }

   public int existeEnlace(String enlazar) {
      for(int i = 0; i < this.enlaces.size(); ++i) {
         Enlace miEnlace = (Enlace)this.enlaces.get(i);
         if (miEnlace.getDestino().equals(enlazar)) {
            return i;
         }
      }

      return -1;
   }

   public double enlacePosicion(int posi) {
      Enlace miEnlace = (Enlace)this.enlaces.get(posi);
      return miEnlace.getPeso();
   }

   public String NodoPosicion(int posi) {
      Enlace miEnlace = (Enlace)this.enlaces.get(posi);
      return miEnlace.getDestino();
   }

   boolean eliminarEnlace(int posicion) {
      if (posicion >= 0 && posicion <= this.enlaces.size()) {
         this.enlaces.remove(posicion);
         --this.enlacesExistentes;
         return true;
      } else {
         JOptionPane.showMessageDialog((Component)null, "No hay enlace en la posicion " + posicion);
         return false;
      }
   }
}