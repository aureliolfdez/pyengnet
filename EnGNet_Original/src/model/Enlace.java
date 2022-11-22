package model;

public class Enlace {
    private String destino;
    private double peso;
 
    public Enlace(String desti, double peso1) {
       this.destino = desti;
       this.peso = peso1;
    }
 
    public Enlace(String desti) {
       this.destino = desti;
       this.peso = -1.0D;
    }
 
    public void modificar(double peso1) {
       this.peso = peso1;
    }
 
    public String getDestino() {
       return this.destino;
    }
 
    public double getPeso() {
       return this.peso;
    }
 }