package model.io;

import java.util.*;

public class Gen {
   private String nombre;
   private ArrayList<Float> experimentos;

   public Gen(String nombre) {
      this.nombre = nombre;
      this.experimentos = new ArrayList();
   }

   public Gen(Gen g) {
      this.setExperimentos(g.getExperimentos());
      this.setNombre(g.getNombre());
   }

   public String getNombre() {
      return this.nombre;
   }

   public void setNombre(String nombre) {
      this.nombre = nombre;
   }

   public ArrayList<Float> getExperimentos() {
      return this.experimentos;
   }

   public void setExperimentos(ArrayList<Float> experimentos) {
      this.experimentos = experimentos;
   }

   public void addExperimento(Float a) {
      this.getExperimentos().add(a);
   }
}