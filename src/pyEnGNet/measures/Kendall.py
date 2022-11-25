import scipy.stats as scp

class Kendall:

    @staticmethod
    def process(dataset, arr1, arr2):
        """
        Recibe dos vectores y calcula el coeficiente de correlación de Kendall entre ambos
            ans -> 1 si pasa el umbral, 0 en caso contrario
            corr -> El valor resultante, en valor absoluto, del coeficiente de correlación de kendall
        :param arr1:
        :param arr2:
        :return:
        """
        ans = 0
        corr, pv = scp.kendalltau(arr1, arr2)
        corr = abs(corr)
        if corr >= dataset.kendall_threshold:
            ans = 1
        return ans, corr    