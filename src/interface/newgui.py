import os
import tkinter as tk
from tkinter import ttk
import customtkinter as ctk

from .checkboxtreeview import CheckboxTreeview

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("dark-blue")


checkboxes = {
    "-CDS-": "CDS",
    "-CENTROMERE-": "centromere",
    "-INTRON-": "intron",
    "-MOBILE_ELEMENT-": "mobile_element",
    "-NCRNA-": "ncRNA",
    "-RRNA-": "rRNA",
    "-TELOMERE-": "telomere",
    "-TRNA-": "tRNA",
    "-3UTR-": "3'UTR",
    "-5UTR-": "5'UTR",
}

__all__ = ["App"]


class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        # treeview Customisation (theme colors are selected)
        bg_color = self._apply_appearance_mode(
            ctk.ThemeManager.theme["CTkFrame"]["fg_color"]
        )
        text_color = self._apply_appearance_mode(
            ctk.ThemeManager.theme["CTkLabel"]["text_color"]
        )
        selected_color = self._apply_appearance_mode(
            ctk.ThemeManager.theme["CTkButton"]["fg_color"]
        )

        treestyle = ttk.Style()
        treestyle.theme_use("default")
        treestyle.configure(
            "Treeview",
            background=bg_color,
            foreground=text_color,
            fieldbackground=bg_color,
            borderwidth=0,
        )
        treestyle.map(
            "Treeview",
            background=[("selected", bg_color)],
            foreground=[("selected", selected_color)],
        )
        self.bind("<<TreeviewSelect>>", lambda event: self.focus_set())

        # configure window
        self.title("BioInfo-Genome")
        self.geometry(f"{1100}x{580}")

        # Frame pour l'arbre de fichiers à gauche
        self.file_tree_frame = ctk.CTkFrame(self, width=450)
        self.file_tree_frame.pack(side="left", fill="y")

        # Frames empilées sur la droite
        self.right_frame = ctk.CTkFrame(self, width=250)
        self.right_frame.pack(side="right", fill="both", expand=False)

        self.checkbox_frame = ctk.CTkFrame(self.right_frame)
        self.checkbox_frame.pack(side="top", fill="both", expand=False)

        # create the checkboxes
        # max 2 checkboxes per line
        self.checkbox_list = []
        self.checkbox_frame_list = []
        for i, (key, value) in enumerate(checkboxes.items()):
            if i % 2 == 0:
                self.checkbox_frame_list.append(ctk.CTkFrame(self.checkbox_frame))
                self.checkbox_frame_list[-1].grid(
                    row=i // 2, column=0, sticky="nsew", padx=5, pady=5
                )
            self.checkbox_list.append(
                ctk.CTkCheckBox(
                    self.checkbox_frame_list[-1],
                    text=value,
                    variable=tk.BooleanVar(),
                    onvalue=True,
                    offvalue=False,
                )
            )
            self.checkbox_list[-1].pack(side="left")

        # Frame pour la prévisualisation au milieu
        self.preview_frame = ctk.CTkFrame(self)
        self.preview_frame.pack(side="top", fill="both", expand=True)

        # Frame tout en bas avec des boutons et du texte
        self.bottom_frame = ctk.CTkFrame(self)
        self.bottom_frame.pack(side="bottom", fill="x")

        # Label pour la prévisualisation
        self.preview_label = ctk.CTkLabel(self.preview_frame, text="Aucune sélection")
        self.preview_label.pack(pady=10)

        # Fonction pour mettre à jour la prévisualisation en fonction de la sélection
        def update_preview(selection):
            # Mettez ici la logique pour mettre à jour la prévisualisation en fonction de la sélection
            self.preview_label.configure(text=f"Prévisualisation pour {selection}")

        # Fonction de gestion de la sélection dans l'arbre de fichiers
        def on_tree_select(event):
            selected_item = self.tree_view.selection()[0]
            update_preview(selected_item)

        # Arbre de fichiers
        self.tree_view = CheckboxTreeview(self.file_tree_frame)
        self.tree_view.pack(expand=True, fill="both")
        self.tree_view["columns"] = ("Depth", "Type", "Size")
        self.tree_view.column("#0", width=300, minwidth=300, stretch=tk.YES)
        self.tree_view.column("Depth", anchor=tk.W, width=50)
        self.tree_view.column("Type", anchor=tk.W, width=50)
        self.tree_view.column("Size", anchor=tk.W, width=50)
        self.tree_view.heading("#0", text="File/Folder")
        self.tree_view.heading("Depth", text="Depth")
        self.tree_view.heading("Type", text="Type")
        self.tree_view.heading("Size", text="Size")
        self.tree_view.bind("<<TreeviewSelect>>", on_tree_select)
        self.fill_tree_view("Results", "")

        # Boutons dans la frame du bas
        self.button1 = ctk.CTkButton(self.bottom_frame, text="Button 1")
        self.button1.grid(row=0, column=0, padx=5, pady=5)

        self.button2 = ctk.CTkButton(self.bottom_frame, text="Button 2")
        self.button2.grid(row=1, column=0, padx=5, pady=5)

        self.button3 = ctk.CTkButton(self.bottom_frame, text="Button 3")
        self.button3.grid(row=1, column=1, padx=5, pady=5)

        self.button4 = ctk.CTkButton(self.bottom_frame, text="Button 4")
        self.button4.grid(row=1, column=2, padx=5, pady=5)

        self.made_with_love_label = ctk.CTkLabel(self.bottom_frame, text="Made with ❤️")
        self.made_with_love_label.grid(row=2, columnspan=2, pady=5)

    def fill_tree_view(self, folder, parent):
        full_path = os.path.join(parent, folder) if parent else folder
        depth = full_path.count(os.sep) if os.sep in full_path else 0
        type_of = "Folder" if os.path.isdir(full_path) else "File"
        size = os.path.getsize(full_path) if os.path.isfile(full_path) else ""

        # insert the current folder or file into the tree view
        self.tree_view.insert(
            parent, "end", full_path, text=folder, values=(depth, type_of, size)
        )

        if os.path.isdir(full_path):
            for item in os.listdir(full_path):
                self.fill_tree_view(item, full_path)


if __name__ == "__main__":
    app = App()
    app.mainloop()
