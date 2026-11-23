Contém os notebooks de química, organizados por tema:

- acidos/01_acidos.ipynb — classificação, ionização e pH de ácidos.
- bases/01_bases.ipynb — hidróxidos, dissociação, pOH e pH.
- sais/01_sais.ipynb — nomenclatura, neutralização, rede iônica e concentração.
- oxidos/01_oxidos.ipynb — classificação e estequiometria da formação de CO2.

Cada notebook é autocontido e inclui teoria, exemplo resolvido, cálculos em
Python, modelo 3D colorido e interativo (arraste para girar 360°), representação
2D, gráfico de dados e exercício com resposta. Compostos iônicos são tratados
como redes/íons, sem inventar moléculas ou ligações covalentes entre os íons.

Abra no Jupyter ou Google Colab e execute as células em ordem. A primeira
célula instala dependências ausentes: rdkit, py3Dmol e matplotlib. Localmente,
use um ambiente virtual com Jupyter/ipykernel. As saídas não são salvas para
manter os arquivos leves; execute as células para gerar as figuras.

O modelo 3D carrega 3Dmol.js de cdn.jsdelivr.net no navegador; requer internet
e confiança no notebook. Se o CDN estiver bloqueado, a representação 2D
continua disponível. Os materiais são teóricos, não roteiros experimentais.
