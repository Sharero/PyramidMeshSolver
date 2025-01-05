# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# # Функция для чтения узлов из файла
# def read_nodes(filename):
#     with open(filename, 'r') as f:
#         nodes = [list(map(float, line.split())) for line in f]
#     return nodes

# # Функция для чтения элементов из файла
# def read_elements(filename):
#     with open(filename, 'r') as f:
#         elements = [list(map(int, line.split())) for line in f]
#     return elements

# # Функция для рисования пирамиды
# def plot_mesh(nodes, elements):
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')

#     # Рисуем элементы
#     for element in elements:
#         # Получаем координаты узлов элемента
#         vertices = [nodes[idx] for idx in element]

#         # Определяем грани пирамиды
#         base = vertices[:4]  # Основание
#         tip = vertices[4]    # Верхушка

#         # Добавляем боковые грани
#         side_face = [base[0], base[1], tip]
#         ax.add_collection3d(Poly3DCollection([side_face], alpha=0.2, edgecolor='k'))
#         side_face = [base[1], base[3], tip]
#         ax.add_collection3d(Poly3DCollection([side_face], alpha=0.2, edgecolor='k'))
#         side_face = [base[2], base[3], tip]
#         ax.add_collection3d(Poly3DCollection([side_face], alpha=0.1, edgecolor='k'))
#         side_face = [base[2], base[0], tip]
#         ax.add_collection3d(Poly3DCollection([side_face], alpha=0.1, edgecolor='k'))
        
#     # Настройка осей
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     ax.autoscale_view()

#     plt.show()

# # Чтение данных из файлов
# nodes = read_nodes('data/nodes.txt')
# elements = read_elements('data/elements.txt')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Define nodes and finite elements for the cube with 6 pyramids
nodes = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
    [0.5, 0.5, 0.5],
]

finite_elements = [
    [0, 4, 1, 5, 8],
    [0, 2, 4, 6, 8],
    [1, 5, 3, 7, 8],
    [2, 3, 6, 7, 8],
    [0, 1, 2, 3, 8],
    [4, 6, 5, 7, 8],
]

# Plotting the cube with pyramids
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Add base pyramids with transparency
for element in finite_elements:
    pyramid = [nodes[element[i]] for i in range(5)]  # Get vertices
    faces = [
        [pyramid[0], pyramid[1], pyramid[4]],
        [pyramid[1], pyramid[3], pyramid[4]],
        [pyramid[3], pyramid[2], pyramid[4]],
        [pyramid[2], pyramid[0], pyramid[4]],
        [pyramid[0], pyramid[1], pyramid[3], pyramid[2]],
    ]
    poly3d = Poly3DCollection(faces, alpha=0.4, linewidths=0.5, edgecolors='k')
    poly3d.set_facecolor("cyan")
    ax.add_collection3d(poly3d)

# Plot node numbers
    ax.text(0,0,0, f'{0}', color='black', fontsize=20)
    ax.text(1,0,0, f'{1}', color='black', fontsize=20)
    ax.text(0,1,0, f'{2}', color='black', fontsize=20)
    ax.text(1,1,0, f'{3}', color='black', fontsize=20)
    ax.text(0,0,1, f'{4}', color='black', fontsize=20)
    ax.text(1,0,1, f'{5}', color='black', fontsize=20)
    ax.text(0,1,1, f'{6}', color='black', fontsize=20)
    ax.text(1,1,1, f'{7}', color='black', fontsize=20)
    ax.text(0.5,0.5,0.5, f'{8}', color='black', fontsize=25)


# Set plot labels and limits
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.set_zlim([0, 1])
ax.view_init(30, 30)

plt.show()
