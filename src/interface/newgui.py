import os
import tkinter as tk
from tkinter import ttk
import customtkinter as ctk

from checkboxtreeview import CheckboxTreeview

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

        # geometry : 2 part left and right and a bottom bar with buttons
        # left part is a treeview
        # right part is a frame with checkboxes and a log window
        # bottom part is a frame with buttons
        self.grid_columnconfigure(0)
        self.grid_columnconfigure(1)
        self.grid_rowconfigure(3, weight=1)

        self.sidebar_frame = ctk.CTkFrame(self, width=300, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, rowspan=4, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(4, weight=1)

        # create the tree view on the left
        self.tree_view = CheckboxTreeview(self.sidebar_frame)
        self.tree_view.pack(fill=tk.BOTH, expand=1)

        # set up columns
        self.tree_view["columns"] = ("Depth", "Type", "Size")
        self.tree_view.column("#0", width=150, minwidth=150, stretch=tk.YES)
        self.tree_view.column("Depth", anchor=tk.W, width=50)
        self.tree_view.column("Type", anchor=tk.W, width=50)
        self.tree_view.column("Size", anchor=tk.W, width=50)

        # add headings
        self.tree_view.heading("#0", text="File/Folder")
        self.tree_view.heading("Depth", text="Depth")
        self.tree_view.heading("Type", text="Type")
        self.tree_view.heading("Size", text="Size")

        # fill tree view with the depth of elements in "Result" folder
        self.fill_tree_view("Results", "")

        # create on the right a frame with
        # - checkboxes on top to select CDS (on multiple rows)
        # - a log window on the bottom that redirects stdout and stderr
        self.checkbox_frame = ctk.CTkFrame(self, corner_radius=0)
        self.checkbox_frame.grid(row=0, column=1, sticky="nsew")
        self.checkbox_frame.grid_columnconfigure(0, weight=1)
        self.checkbox_frame.grid_rowconfigure(0, weight=1)

        # create the checkboxes
        # max 4 checkboxes per line
        self.checkboxes = {}
        self.checkbox_rows = []
        self.checkbox_columns = []
        self.checkbox_frame.grid_columnconfigure(0, weight=1)
        self.checkbox_frame.grid_rowconfigure(0, weight=1)

        for i, (key, value) in enumerate(checkboxes.items()):
            self.checkboxes[key] = ctk.CTkCheckBox(self.checkbox_frame, text=value)
            self.checkboxes[key].grid(row=i // 4, column=i % 4, sticky="w")
            self.checkbox_rows.append(i // 4)
            self.checkbox_columns.append(i % 4)

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
