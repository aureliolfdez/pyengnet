package util;

import java.io.*;
import java.util.*;

public class ReadFromPropertiesSingleton {
   private static ReadFromPropertiesSingleton instance = null;
   private Properties prop = new Properties();

   private ReadFromPropertiesSingleton() {
      FileInputStream input = null;

      try {
         input = new FileInputStream("config.properties");
         this.prop.load(input);
         System.out.println("Config file loaded");
      } catch (IOException var11) {
         var11.printStackTrace();
      } finally {
         if (input != null) {
            try {
               input.close();
            } catch (IOException var10) {
               var10.printStackTrace();
            }
         }

      }

   }

   public static ReadFromPropertiesSingleton getInstance() {
      if (instance == null) {
         instance = new ReadFromPropertiesSingleton();
      }

      return instance;
   }

   public Properties getProp() {
      return this.prop;
   }
}