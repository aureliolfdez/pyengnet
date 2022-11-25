import scipy.stats as scp

class Spearman:
    
    @staticmethod
    def process(dataset, arr1, arr2):
        """
        Recibe dos vectores y calcula el coeficiente de correlación de spearman entre ellos, devuelve:
            ans -> 1 si pasa el umbral, 0 en caso contrario
            corr -> El coeficiente de correlación de spearman en valor absoluto
        :param arr1:
        :param arr2:
        :return:
        """
        ans = 0
        corr, pv = scp.spearmanr(arr1, arr2)
        corr = abs(corr)
        if corr >= dataset.spearman_threshold:
            ans = 1
        return ans, corr