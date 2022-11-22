package model.performance;

public class Arista {
    private String origen;
    private String destino;
    private float etiqueta;
    private float valorMate;
    private float valorBio;
    private float NonDubiousDegree;
    private double biologicalSimilitude = 0.0D;
    private double biologicalDisimilitude = 0.0D;
 
    public float getNonDubiousDegree() {
       return this.NonDubiousDegree;
    }
 
    public void setNonDubiousDegree(float nonDubiousDegree) {
       this.NonDubiousDegree = nonDubiousDegree;
    }
 
    public double getBiologicalSimilitude() {
       return this.biologicalSimilitude;
    }
 
    public void setBiologicalSimilitude(double biologicalSimilitude) {
       this.biologicalSimilitude = biologicalSimilitude;
    }
 
    public double getBiologicalDisimilitude() {
       return this.biologicalDisimilitude;
    }
 
    public void setBiologicalDisimilitude(double biologicalDisimilitude) {
       this.biologicalDisimilitude = biologicalDisimilitude;
    }
 
    public float getValorMate() {
       return this.valorMate;
    }
 
    public void setValorMate(float valorMate) {
       this.valorMate = valorMate;
    }
 
    public float getValorBio() {
       return this.valorBio;
    }
 
    public void setValorBio(float valorBio) {
       this.valorBio = valorBio;
    }
 
    public Arista(String nombre, String nombre2, float etiqueta) {
       this.origen = nombre;
       this.destino = nombre2;
       this.etiqueta = etiqueta;
    }
 
    public float getEtiqueta() {
       return this.etiqueta;
    }
 
    public void setEtiqueta(float etiqueta) {
       this.etiqueta = etiqueta;
    }
 
    public String getOrigen() {
       return this.origen;
    }
 
    public void setOrigen(String origen) {
       this.origen = origen;
    }
 
    public String getDestino() {
       return this.destino;
    }
 
    public void setDestino(String destino) {
       this.destino = destino;
    }
 
    public boolean equals(Object obj) {
       Arista a = (Arista)obj;
       return (a.getOrigen().equals(this.getOrigen()) || a.getDestino().equals(this.getOrigen())) && (a.getOrigen().equals(this.getDestino()) || a.getDestino().equals(this.getDestino()));
    }
 
    public int hashCode() {
       return 1;
    }
 }
