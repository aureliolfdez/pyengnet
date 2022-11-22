package model;

public class Arco {
    private String inicial;
    private String terminal;
    private float peso;
    private float peso1;
    private float peso2;
    private float peso3;
 
    public Arco(String ini, String ter, float pes) {
       this.inicial = ini;
       this.terminal = ter;
       this.peso = pes;
    }
 
    public Arco(String ini, String ter, float pes1, float pes2, float pes3) {
       this.inicial = ini;
       this.terminal = ter;
       this.peso1 = pes1;
       this.peso2 = pes2;
       this.peso3 = pes3;
       this.peso = (pes1 + pes2 + pes3) / 3.0F;
    }
 
    public float getPeso() {
       return this.peso;
    }
 
    public void setPeso(float peso) {
       this.peso = peso;
    }
 
    public String getInicial() {
       return this.inicial;
    }
 
    public void setInicial(String inicial) {
       this.inicial = inicial;
    }
 
    public String getTerminal() {
       return this.terminal;
    }
 
    public void setTerminal(String terminal) {
       this.terminal = terminal;
    }
 
    public String toString() {
       return this.inicial + "\t" + this.terminal + "\t" + this.peso;
    }
 
    public String toStringESM() {
       return this.inicial + "\t" + this.terminal + "\t" + this.peso1 + "\t" + this.peso2 + "\t" + this.peso3;
    }
 
    public boolean igualesAsociacion(Arco a) {
       return this.getInicial().equals(a.getInicial()) && this.getTerminal().equals(a.getTerminal()) || this.getInicial().equals(a.getTerminal()) && this.getTerminal().equals(a.getInicial());
    }
 
    public boolean equals(Object obj) {
       Arco a = (Arco)obj;
       return this.getInicial().equals(a.getInicial()) && this.getTerminal().equals(a.getTerminal()) || this.getInicial().equals(a.getTerminal()) && this.getTerminal().equals(a.getInicial());
    }
 }